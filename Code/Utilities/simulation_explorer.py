"""

Interactive explorer for simulation results stored in HDF5 (.h5) files.

The window has a menu on the left and the figure on the right. 

In the menu:

1) "Open folder..." selects the folder where the simulations are stored
   (the .h5 files can be directly inside it or spread over subfolders);
   Only chemistry simulations (files with a '/chemistry' group) are considered;
   "Simulations" then shows its subfolders as a tree, all ticked: unticking a folder
   leaves its simulations out of everything below (parameters, plots and subspecies);

2) "Parameter" chooses, among the attributes common to all of them, in function
   of which one the results are represented;

3) "Plot" chooses what to plot (Reduced Field, Final Densities, Populations or
   Particle Balance) and, for densities/populations/balance, which subspecies to plot in the same figure.
   An option is offered if at least one file has its data, and it uses only those files.

The figure is drawn as soon as everything is selected; its toolbar allows to zoom,
edit the title, labels and axis scales, and save it.

For Final Densities, "Animate in time" shows instead how the densities evolve in time,
read from '/chemistry/chemSolutionTime' at a chosen number of instants (evenly spaced
between a start and an end time; each file uses its nearest time step). 
The animation plays in the figure and "Save GIF..." saves it.

The same questions can also be answered in the terminal (e.g. over ssh, without a
display); the figure is then shown in a window or saved to a file.

Run with:

python simulation_explorer.py          (window)
python simulation_explorer.py --cli    (terminal)

"""

__author__ = 'Pedro Pereira <pedro.t.pereira@tecnico.ulisboa.pt>'
__maintainer__ = 'Pedro Pereira <pedro.t.pereira@tecnico.ulisboa.pt>'
__date__ = '13 November 2026'
__version__ = 'Revision: 0.4.0'

import argparse
import os
import re
import shutil
import sys
import textwrap

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib import animation


# Shown when the tool starts: in the terminal before anything is asked, in the window at the
# bottom left, and in --help
INTRO_TEXT = ("This is a tool to plot the outputs of LoKI-GM simulations, representing the results "
              "as a function of the working condition parameters common to all of them. The selected "
              "folder should contain chemistry simulations, with their .h5 output files.")

PLOT_OPTIONS = [
    ('Reduced Field', '/electronKinetics/reducedField'),
    ('Final Densities', '/chemistry/finalDensities'),
    ('Populations', '/chemistry/finalDensities'),
    ('Particle Balance', '/chemistry/particleBalance'),
]

FIELD_PATH = '/electronKinetics/reducedField'
DENSITIES_PATH = '/chemistry/finalDensities'
BALANCE_PATH = '/chemistry/particleBalance'
TIME_PATH = '/chemistry/chemSolutionTime'
CHEMISTRY_PATH = '/chemistry'
EXCLUDED_BALANCE_SUBSPECIES = {'N2(X)'}
SUBSPECIES_PLOTS = ('Final Densities', 'Populations', 'Particle Balance')


def find_h5_files(path):
    # Walks through every folder inside path and collects the .h5 files
    h5_files = []
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if filename.endswith('.h5'):
                h5_files.append(os.path.join(dirpath, filename))
    return sorted(h5_files)


def tree_branches(path, h5_files):
    # For each file, the labels of its branch in the "Simulations" tree, below the opened folder
    # (e.g. ['Test_Run', 'sim_1']). A folder with a single .h5 file and no simulations in its
    # subfolders stands for that file; otherwise the file names are shown under the folder
    folders = [os.path.normpath(os.path.relpath(os.path.dirname(f), path)) for f in h5_files]
    folders = ['' if folder == '.' else folder for folder in folders]
    n_files = {folder: folders.count(folder) for folder in folders}

    def has_subfolders(folder):
        prefix = folder + os.sep if folder else ''
        return any(other != folder and other.startswith(prefix) for other in n_files)

    branches = []
    for file_path, folder in zip(h5_files, folders):
        parts = folder.split(os.sep) if folder else []
        if not parts or n_files[folder] > 1 or has_subfolders(folder):
            parts.append(os.path.basename(file_path))
        branches.append(parts)
    return branches


def to_python_value(value):
    # Converts an HDF5 attribute value into a hashable python value
    if isinstance(value, bytes):
        return value.decode('utf-8')
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return to_python_value(value.reshape(-1)[0])
        return tuple(to_python_value(v) for v in value.reshape(-1))
    if isinstance(value, np.generic):
        return value.item()
    return value


def read_file_info(file_path):
    # Reads the root attributes of an .h5 file and which plot groups it contains;
    # returns None if it is not a chemistry simulation (LoKI-GM creates the /chemistry group
    # at the start of the run only when the chemistry is enabled)
    with h5py.File(file_path, 'r') as f:
        if CHEMISTRY_PATH not in f:
            return None
        attrs = {name: to_python_value(f.attrs[name]) for name in f.attrs}
        groups = {name for name, path in PLOT_OPTIONS if path in f}
    return attrs, groups


def read_folder(path):
    # Reads the info of the chemistry simulations inside path, skipping the files that cannot
    # be opened; also returns how many .h5 files were ignored for having no chemistry
    files, all_attrs, all_groups, n_without_chemistry = [], [], [], 0
    for file_path in find_h5_files(path):
        try:
            info = read_file_info(file_path)
        except OSError as e:
            print(f"Warning: could not read '{file_path}' ({e}). Skipping.")
            continue
        if info is None:
            n_without_chemistry += 1
            continue
        files.append(file_path)
        all_attrs.append(info[0])
        all_groups.append(info[1])
    return files, all_attrs, all_groups, n_without_chemistry


def folder_summary(n_files):
    # e.g. "14 chemistry simulation(s)"
    text = f"{n_files} chemistry simulation(s)"

    return text


def common_attributes(all_attrs):
    # Returns the attributes present in every file (order of the first file)
    # and, for each one, the set of distinct values found
    common_names = set(all_attrs[0])
    for attrs in all_attrs[1:]:
        common_names &= set(attrs)

    ordered = [name for name in all_attrs[0] if name in common_names]
    distinct = {name: set(attrs[name] for attrs in all_attrs) for name in ordered}
    return ordered, distinct


def available_plot_options(all_groups):
    # The plot options whose group exists in at least one file
    return [name for name, _ in PLOT_OPTIONS if any(name in groups for groups in all_groups)]


def files_for_plot(plot, h5_files, all_attrs, all_groups):
    # Keeps only the files (and their attributes) that have the data of the chosen plot,
    # plus a message saying how many were left out ("" if none)
    keep = [i for i, groups in enumerate(all_groups) if plot in groups]
    n_out = len(h5_files) - len(keep)
    message = (f"{n_out} of {len(h5_files)} simulation(s) not used, for having no {plot} data." if n_out else "")
    return [h5_files[i] for i in keep], [all_attrs[i] for i in keep], message


def list_subspecies(h5_files, plot):
    # Returns the (parent, subspecies) pairs found in any file, in order of first appearance
    pairs = []
    seen = set()

    for file_path in h5_files:
        with h5py.File(file_path, 'r') as f:
            densities = f[DENSITIES_PATH]

            if plot == 'Particle Balance':
                balance = f[BALANCE_PATH]
                parents = [p for p in balance.keys() if balance[p].shape[0] > 0]
            else:
                parents = list(densities.keys())

            for parent in parents:
                if parent not in densities:
                    continue
                dataset = densities[parent]
                if dataset.ndim < 2 or dataset.shape[-1] == 0:
                    continue

                n_sub = dataset.shape[1]
                if plot == 'Particle Balance':
                    n_sub = min(n_sub, f[BALANCE_PATH][parent].shape[0])

                for i in range(n_sub):
                    name = dataset[0, i]['Species'].decode('utf-8')
                    if plot == 'Particle Balance' and name in EXCLUDED_BALANCE_SUBSPECIES:
                        continue
                    if (parent, name) not in seen:
                        seen.add((parent, name))
                        pairs.append((parent, name))

    return pairs


def ask_choice(question, labels):
    # Prints a numbered list and asks the user to choose one entry by its number
    print(f"\n{question}")
    for i, label in enumerate(labels, start=1):
        print(f"  {i:>3}) {label}")

    while True:
        choice = input(f"Choice number (1-{len(labels)}): ").strip()
        if choice.isdigit() and 1 <= int(choice) <= len(labels):
            return int(choice) - 1
        print(f"Invalid choice, please enter a number between 1 and {len(labels)}.")


def parse_selection(text, n_options):
    # Parses "1 4 7-9" or "2,5" into a list of 0-based indices; returns None if invalid
    tokens = [t for t in re.split(r'[,\s]+', text.strip()) if t]
    if not tokens:
        return None

    indices = []
    for token in tokens:
        range_match = re.fullmatch(r'(\d+)-(\d+)', token)
        if range_match:
            start, end = int(range_match.group(1)), int(range_match.group(2))
            if start > end:
                return None
            numbers = range(start, end + 1)
        elif token.isdigit():
            numbers = [int(token)]
        else:
            return None

        for number in numbers:
            if not 1 <= number <= n_options:
                return None
            if number - 1 not in indices:
                indices.append(number - 1)

    return indices


def ask_multiple_choice(question, pairs):
    # Prints the (parent, subspecies) pairs grouped by parent and asks for several numbers
    print(f"\n{question}")
    current_parent = None
    for i, (parent, name) in enumerate(pairs, start=1):
        if parent != current_parent:
            print(f"  --- {parent} ---")
            current_parent = parent
        print(f"  {i:>4}) {name}")

    while True:
        text = input(f"Numbers separated by spaces or commas, ranges allowed (e.g. 1 4 7-9): ")
        indices = parse_selection(text, len(pairs))
        if indices is not None:
            return [pairs[i] for i in indices]
        print(f"Invalid selection, please use numbers between 1 and {len(pairs)}.")


def plot_reduced_field(ax, h5_files, all_attrs, parameter):
    # Plots the reduced electric field E/N [Td] in function of the chosen parameter
    # (as in field_analyser of densities_field_cases.py); returns a message if there is nothing to plot
    x_values = []
    y_values = []
    for file_path, attrs in zip(h5_files, all_attrs):
        with h5py.File(file_path, 'r') as f:
            field = f[FIELD_PATH]
            if field.shape[-1] == 0:
                continue
            x_values.append(attrs[parameter])
            y_values.append(float(field[0, 0]))

    if not x_values:
        return "No reduced field information found in the dataset."

    order = np.argsort(x_values)
    ax.plot(np.asarray(x_values)[order], np.asarray(y_values)[order], marker='o')
    ax.set_xlabel(parameter, fontsize=12)
    ax.set_ylabel("E/N (Td)", fontsize=12)
    ax.tick_params(labelsize=12)
    return ""


def figure_filename(kind, parameter, extra=()):
    # Builds e.g. "sim_field_Discharge_current.png": units between parenthesis are dropped
    # from the parameter, and any character not valid in a file name becomes "_"
    parameter_tag = re.sub(r'\([^)]*\)', '', parameter)
    parts = [re.sub(r'[^\w\-.]+', '_', part).strip('_')
             for part in ("sim", kind, parameter_tag, *extra)]
    return "_".join(part for part in parts if part) + ".png"


def read_density(f, parent, name, column='Abs.Density'):
    # Absolute density ('Abs.Density') or population ('Population') of a subspecies,
    # found by name (None if not in this file)
    if parent not in f[DENSITIES_PATH]:
        return None
    dataset = f[DENSITIES_PATH][parent]
    if dataset.ndim < 2 or dataset.shape[-1] == 0:
        return None
    for i in range(dataset.shape[1]):
        row = dataset[0, i]
        if row['Species'].decode('utf-8') == name:
            return float(row[column])
    return None


def density_ylabel(subspecies):
    # Wall species are surface densities (m-2), the others volume densities (m-3)
    units = {'m$^{-2}$' if parent.lower().startswith('wall') else 'm$^{-3}$'
             for parent, _ in subspecies}
    return f"Absolute density ({' / '.join(sorted(units))})"


def plot_final_densities(ax, h5_files, all_attrs, parameter, subspecies, populations=False):
    # Plots the absolute density (or, with populations=True, the population relative to the
    # parent species) of several subspecies in function of the chosen parameter, all in the
    # same figure (as densities_analyser with -hold in densities_field_cases.py)
    column, quantity = ('Population', 'population') if populations else ('Abs.Density', 'density')
    data = {pair: ([], []) for pair in subspecies}
    for file_path, attrs in zip(h5_files, all_attrs):
        with h5py.File(file_path, 'r') as f:
            for parent, name in subspecies:
                value = read_density(f, parent, name, column)
                if value is None:
                    continue
                data[(parent, name)][0].append(attrs[parameter])
                data[(parent, name)][1].append(value)

    if not any(x for x, _ in data.values()):
        return f"No {quantity} information found in the dataset for the selected subspecies."

    missing = []
    for (parent, name), (x_values, y_values) in data.items():
        if not x_values:
            missing.append(name)
            continue
        order = np.argsort(x_values)
        ax.plot(np.asarray(x_values)[order], np.asarray(y_values)[order], marker='o', label=name)

    ax.set_xlabel(parameter, fontsize=12)
    ax.set_ylabel("Population" if populations else density_ylabel(subspecies), fontsize=12)
    ax.tick_params(labelsize=12)
    ax.yaxis.get_offset_text().set_fontsize(12)  # the "1e22" multiplier
    ax.legend(fontsize=12)
    return f"No {quantity} information found for {', '.join(missing)}." if missing else ""


def time_range(h5_files):
    # Returns the smallest and the largest time in chemSolutionTime over all the files
    # (None if no file has it), to suggest the limits of the animation
    t_min, t_max = np.inf, -np.inf
    for file_path in h5_files:
        with h5py.File(file_path, 'r') as f:
            if TIME_PATH not in f or f[TIME_PATH].shape[-1] == 0:
                continue
            times = f[TIME_PATH].fields('Time')[0]
            t_min, t_max = min(t_min, times.min()), max(t_max, times.max())
    return (t_min, t_max) if t_max > t_min else None


def time_instants(t_start, t_end, n_frames):
    # Evenly spaced time instants of the animation; returns (instants, "") or (None, error message)
    if not t_end > t_start:
        return None, "The end time must be larger than the start time."
    if n_frames < 2:
        return None, "The animation needs at least 2 frames."
    return np.linspace(t_start, t_end, n_frames), ""


def read_time_densities(h5_files, all_attrs, parameter, subspecies, instants):
    # Reads the densities of the subspecies from chemSolutionTime at each instant (each file
    # uses its time step nearest to the instant); every file is opened only once.
    # Returns frames[k][(parent, name)] = (x values, y values), sorted by x
    frames = [{pair: ([], []) for pair in subspecies} for _ in instants]
    for file_path, attrs in zip(h5_files, all_attrs):
        with h5py.File(file_path, 'r') as f:
            if TIME_PATH not in f or f[TIME_PATH].shape[-1] == 0:
                continue
            dataset = f[TIME_PATH]
            pairs = [pair for pair in subspecies if pair[1] in dataset.dtype.names]
            if not pairs:
                continue
            table = dataset.fields(['Time'] + list(dict.fromkeys(name for _, name in pairs)))[0]
            steps = np.abs(table['Time'][None, :] - np.asarray(instants)[:, None]).argmin(axis=1)
            for frame, step in zip(frames, steps):
                for pair in pairs:
                    frame[pair][0].append(attrs[parameter])
                    frame[pair][1].append(float(table[step][pair[1]]))

    for frame in frames:
        for pair, (x_values, y_values) in frame.items():
            order = np.argsort(x_values)
            frame[pair] = (np.asarray(x_values)[order], np.asarray(y_values)[order])
    return frames


def make_density_animation(fig, ax, frames, instants, parameter, subspecies):
    # Animates the densities of several subspecies in function of the chosen parameter, one
    # frame per time instant (as densities_cases_time_animation.py). The axes limits cover
    # all the frames, so they don't jump. Returns (animation, message); the animation must
    # be kept in a variable, otherwise it is deleted and stops.
    all_x = np.concatenate([x for frame in frames for x, _ in frame.values()])
    all_y = np.concatenate([y for frame in frames for _, y in frame.values()])
    if all_x.size == 0:
        return None, "No density information found in chemSolutionTime for the selected subspecies."

    missing = [name for parent, name in subspecies
               if not any(frame[(parent, name)][0].size for frame in frames)]
    lines = {pair: ax.plot([], [], marker='o', label=pair[1])[0] for pair in subspecies}
    ax.update_datalim(np.column_stack([all_x, all_y]))
    ax.autoscale_view()
    ax.set_xlabel(parameter, fontsize=12)
    ax.set_ylabel(density_ylabel(subspecies), fontsize=12)
    ax.tick_params(labelsize=12)
    ax.yaxis.get_offset_text().set_fontsize(12)  # the "1e22" multiplier
    ax.legend(fontsize=12)

    def update(k):
        for pair, line in lines.items():
            line.set_data(*frames[k][pair])
        ax.set_title(f"t = {instants[k]:.3g} s", fontsize=12)
        return list(lines.values())

    update(0)
    anim = animation.FuncAnimation(fig, update, frames=len(instants), interval=500, repeat=True)
    message = f"No density information found in chemSolutionTime for {', '.join(missing)}." if missing else ""
    return anim, message


def save_animation(anim, path, fps=2):
    # Saves the animation as a GIF; the 'tight' bounding box used for the figures would give
    # frames of different sizes, so the standard one is used here
    with plt.rc_context({'savefig.bbox': 'standard'}):
        anim.save(path, writer='pillow', fps=fps, dpi=100)


def ask_number(question, minimum=None, integer=False):
    # Asks for a number (re-asks if invalid or below minimum)
    while True:
        text = input(question).strip()
        try:
            value = int(text) if integer else float(text)
        except ValueError:
            print("Please enter " + ("an integer." if integer else "a number."))
            continue
        if minimum is not None and value < minimum:
            print(f"Please enter a number not smaller than {minimum}.")
            continue
        return value


def ask_limit():
    # Asks for the threshold (%) from which reactions are shown in the particle balance
    while True:
        text = input("\nThreshold (%) from which reactions are considered (0-100): ").strip()
        try:
            limit = float(text)
        except ValueError:
            print("Please enter a number.")
            continue
        if 0 <= limit <= 100:
            return limit
        print("Please specify a limit between 0 and 100.")


def read_balance(f, parent, name):
    # Returns the (reaction, contribution %) pairs for the gains (positive values) and
    # losses (negative values) of a subspecies, or None if not in this file
    if parent not in f[BALANCE_PATH] or parent not in f[DENSITIES_PATH]:
        return None
    balance = f[BALANCE_PATH][parent]
    if balance.shape[0] == 0 or balance.shape[-1] == 0:
        return None

    # The particleBalance dataset follows the order of finalDensities, without N2(X)
    dataset = f[DENSITIES_PATH][parent]
    names = [dataset[0, i]['Species'].decode('utf-8') for i in range(dataset.shape[1])]
    names = [n for n in names if n not in EXCLUDED_BALANCE_SUBSPECIES]
    if name not in names or names.index(name) >= balance.shape[0]:
        return None

    child = balance[names.index(name), :, 0, :]
    entries = []
    for row, sign, kind in ((child[0, :], 1.0, 'gain'), (child[1, :], -1.0, 'loss')):
        for reaction, value in zip(row['Reactions'], row['Contribution']):
            entries.append((reaction.decode('utf-8'), kind, sign * float(value)))
    return entries


def plot_particle_balance(ax, h5_files, all_attrs, parameter, subspecies, limit):
    # Plots the contribution (%) of the main creation (positive) and destruction (negative)
    # mechanisms of the selected subspecies, in function of the chosen parameter
    # (as in balance_drive_conditions.py: if a reaction is above the limit at any point,
    # it is plotted for all the data points)
    multiple = len(subspecies) > 1
    all_data = {}  # (subspecies, reaction, gain/loss) -> ([x], [y]), for every reaction
    for file_path, attrs in zip(h5_files, all_attrs):
        with h5py.File(file_path, 'r') as f:
            for parent, name in subspecies:
                entries = read_balance(f, parent, name)
                if entries is None:
                    continue
                for reaction, kind, value in entries:
                    x_list, y_list = all_data.setdefault((name, reaction, kind), ([], []))
                    x_list.append(attrs[parameter])
                    y_list.append(value)

    # Keep only the reactions that are above the limit in at least one point
    data = {key: points for key, points in all_data.items()
            if any(abs(value) > limit for value in points[1])}

    if not data:
        return f"No reactions above {limit}% found for the selected subspecies."

    # The default palette has 10 colours; with more lines use tab20 so they don't repeat
    n_default = len(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    colors = plt.cm.tab20.colors if len(data) > n_default else None

    for i, ((name, reaction, kind), (x_values, y_values)) in enumerate(data.items()):
        order = np.argsort(x_values)
        label = f"{name}: {reaction}" if multiple else reaction
        color = colors[i % len(colors)] if colors else None
        ax.plot(np.asarray(x_values)[order], np.asarray(y_values)[order], marker='o', label=label,
                color=color)

    ylabel = "Reaction contributions (%)" if multiple else f"{subspecies[0][1]} reaction contributions (%)"
    ax.set_xlabel(parameter, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.tick_params(labelsize=12)
    # Placed automatically inside the plot, default font size; it can be moved with the mouse
    ax.legend().set_draggable(True)
    return ""


def gui():
    # Window with the menu on the left and the figure on the right; matplotlib's toolbar,
    # above the figure, allows to zoom, edit the labels and axis scales, and save it.
    # Qt is imported here so that the command line also works where Qt is not available
    from matplotlib.backends.qt_compat import QtCore, QtGui, QtWidgets
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
    from matplotlib.figure import Figure

    app = QtWidgets.QApplication(sys.argv)

    window = QtWidgets.QWidget(windowTitle="Simulation explorer")
    folder_button = QtWidgets.QPushButton("Open folder...")
    folder_label = QtWidgets.QLabel("No folder selected")
    # Subfolders of the opened folder, with a tick box each: only the ticked simulations are used
    folder_tree = QtWidgets.QTreeWidget(headerHidden=True)
    folder_tree.setMaximumHeight(180)  # so the subspecies list keeps most of the space
    all_button = QtWidgets.QPushButton("All")
    none_button = QtWidgets.QPushButton("None")
    tree_widget = QtWidgets.QWidget()
    tree_layout = QtWidgets.QGridLayout(tree_widget)
    tree_layout.setContentsMargins(0, 0, 0, 0)
    tree_layout.addWidget(folder_tree, 0, 0, 1, 2)
    tree_layout.addWidget(all_button, 1, 0)
    tree_layout.addWidget(none_button, 1, 1)
    parameter_box = QtWidgets.QComboBox()  # long names are cut with "..." but shown complete when opened
    parameter_box.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
    plot_box = QtWidgets.QComboBox()
    subspecies_list = QtWidgets.QListWidget()  # several can be selected, one click each
    subspecies_list.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)
    limit_box = QtWidgets.QDoubleSpinBox(maximum=100, value=5, suffix=' %')
    animate_check = QtWidgets.QCheckBox("Animate in time")
    # Time settings of the animation, grouped so they can be shown or hidden together
    time_validator = QtGui.QDoubleValidator(bottom=0)
    time_validator.setLocale(QtCore.QLocale.c())  # "." as decimal separator, as float() expects
    t_start_box = QtWidgets.QLineEdit(validator=time_validator)
    t_end_box = QtWidgets.QLineEdit(validator=time_validator)
    n_frames_box = QtWidgets.QSpinBox(minimum=2, maximum=500, value=20)
    n_frames_box.setKeyboardTracking(False)  # redraws when the number is complete, not per digit
    time_widget = QtWidgets.QWidget()
    time_form = QtWidgets.QFormLayout(time_widget)
    time_form.setContentsMargins(0, 0, 0, 0)
    time_form.addRow("Start time (s)", t_start_box)
    time_form.addRow("End time (s)", t_end_box)
    time_form.addRow("Number of frames", n_frames_box)
    save_gif_button = QtWidgets.QPushButton("Save GIF...", enabled=False)
    clear_button = QtWidgets.QPushButton("Clear plot")
    status = QtWidgets.QLabel(wordWrap=True)

    menu_widget = QtWidgets.QWidget()
    menu_widget.setFixedWidth(300)  # so it does not change when rows are shown or hidden
    menu = QtWidgets.QFormLayout()
    menu.setContentsMargins(0, 0, 0, 0)
    menu.setRowWrapPolicy(QtWidgets.QFormLayout.RowWrapPolicy.WrapAllRows)  # labels above the widgets
    menu.addRow(folder_button)
    menu.addRow(folder_label)
    menu.addRow("Simulations", tree_widget)
    menu.addRow("Parameter (x axis)", parameter_box)
    menu.addRow("Plot (y axis)", plot_box)
    menu.addRow("Subspecies", subspecies_list)
    menu.addRow("Show reactions above", limit_box)
    menu.addRow(animate_check)
    menu.addRow(time_widget)
    menu.addRow(save_gif_button)
    menu.addRow(clear_button)
    menu.addRow(status)

    # The introduction stays at the bottom of the menu; the menu above takes the free space
    intro_label = QtWidgets.QLabel(INTRO_TEXT, wordWrap=True)
    intro_label.setForegroundRole(QtGui.QPalette.ColorRole.PlaceholderText)  # grey, in light and dark themes
    menu_column = QtWidgets.QVBoxLayout(menu_widget)
    menu_column.setContentsMargins(0, 0, 0, 0)
    menu_column.addLayout(menu, stretch=1)
    menu_column.addWidget(intro_label)

    figure = Figure()
    canvas = FigureCanvasQTAgg(figure)
    toolbar = NavigationToolbar2QT(canvas, window)
    plot_area = QtWidgets.QVBoxLayout()
    plot_area.addWidget(toolbar)
    plot_area.addWidget(canvas)

    layout = QtWidgets.QHBoxLayout(window)
    layout.addWidget(menu_widget)
    layout.addLayout(plot_area, stretch=1)

    # all_*: every chemistry simulation of the opened folder; files/attrs/groups: the ticked ones
    data = {'all_files': [], 'all_attrs': [], 'all_groups': [], 'files': [], 'attrs': [], 'groups': [],
            'pairs': [], 'path': '', 'anim': None}

    def plot_files(plot):
        # Files with the data of the chosen plot, their attributes and a note on the others
        return files_for_plot(plot, data['files'], data['attrs'], data['groups'])

    def show_options(loaded):
        # The options are shown only once a folder is read and there are files to plot
        for widget in (parameter_box, plot_box, clear_button):
            menu.setRowVisible(widget, loaded)

    def open_folder():
        # Reads the chosen folder, lists its simulations in the tree (all ticked) and fills the menus
        path = QtWidgets.QFileDialog.getExistingDirectory(window, "Select the folder with the simulations")
        if not path:  # cancelled
            return
        files, all_attrs, all_groups, _ = read_folder(path)
        if not files:
            # The folder that was open (if any) stays, so the message says which one is shown
            message = f"No chemistry simulations with .h5 files found in '{path}'."
            if data['all_files']:
                message += f" Still showing '{os.path.basename(os.path.normpath(data['path']))}'."
            status.setText(message)
            return
        data.update(all_files=files, all_attrs=all_attrs, all_groups=all_groups, path=path)
        plt.rcParams['savefig.directory'] = path  # where the save dialog starts
        build_tree(path, files)
        menu.setRowVisible(tree_widget, True)

        # Suggests the whole simulated time for the animation
        times = time_range(files)
        animate_check.setEnabled(times is not None)
        animate_check.setToolTip("" if times else "The files have no chemSolutionTime data.")
        if times is None:
            animate_check.setChecked(False)
        else:
            t_start_box.setText(f"{times[0]:.3g}")
            t_end_box.setText(f"{times[1]:.3g}")

        apply_selection()

    def build_tree(path, files):
        # One item per folder; the items of the simulations keep the index of their file
        folder_tree.blockSignals(True)  # so building it does not update the plot at each item
        folder_tree.clear()
        root = QtWidgets.QTreeWidgetItem(folder_tree, [os.path.basename(os.path.normpath(path))])
        items = {(): root}
        for index, branch in enumerate(tree_branches(path, files)):
            for depth in range(1, len(branch) + 1):
                key = tuple(branch[:depth])
                if key not in items:
                    items[key] = QtWidgets.QTreeWidgetItem(items[key[:-1]], [branch[depth - 1]])
            items[tuple(branch)].setData(0, QtCore.Qt.ItemDataRole.UserRole, index)
            items[tuple(branch)].setToolTip(0, files[index])
        for item in items.values():
            # A folder is ticked when all inside it are, and ticking it ticks (or unticks) them all
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable
                          | QtCore.Qt.ItemFlag.ItemIsAutoTristate)
            item.setCheckState(0, QtCore.Qt.CheckState.Checked)
        folder_tree.expandAll()
        folder_tree.blockSignals(False)

    def ticked_indices():
        # Indices (in data['all_files']) of the ticked simulations
        indices = []
        iterator = QtWidgets.QTreeWidgetItemIterator(folder_tree)
        while iterator.value():
            item = iterator.value()
            index = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
            if index is not None and item.checkState(0) == QtCore.Qt.CheckState.Checked:
                indices.append(index)
            iterator += 1
        return sorted(indices)

    def set_all_ticked(ticked):
        # The "All" and "None" buttons: ticking the top folder ticks everything inside it
        root = folder_tree.topLevelItem(0)
        if root is not None:
            root.setCheckState(0, QtCore.Qt.CheckState.Checked if ticked else QtCore.Qt.CheckState.Unchecked)

    # Ticking a folder changes the state of every item inside it, one signal each: the selection
    # is applied once, after all of them (a timer of 0 ms runs when Qt has finished with the click)
    selection_timer = QtCore.QTimer(singleShot=True, interval=0)

    def apply_selection():
        # Fills the menus with the parameters and plot options of the ticked simulations,
        # keeping the chosen parameter, plot and subspecies when they are still available
        indices = ticked_indices()
        data.update(files=[data['all_files'][i] for i in indices],
                    attrs=[data['all_attrs'][i] for i in indices],
                    groups=[data['all_groups'][i] for i in indices])
        summary = folder_summary(len(indices))
        if len(indices) < len(data['all_files']):
            summary += f" selected, of {len(data['all_files'])}"
        folder_label.setText(f"{os.path.basename(os.path.normpath(data['path']))}: {summary}")

        # The last parameter and plot chosen (also after "None", when the menus were empty)
        previous_parameter = parameter_box.currentData() or data.get('last_parameter')
        previous_plot = plot_box.currentData() or data.get('last_plot')
        data.update(last_parameter=previous_parameter, last_plot=previous_plot)
        for box in (parameter_box, plot_box):  # refilled without updating the plot at each step
            box.blockSignals(True)
            box.clear()

        if indices:
            names, distinct = common_attributes(data['attrs'])
            for name in names:
                n_values = len(distinct[name])
                tag = f"varies: {n_values} values" if n_values > 1 else f"constant = {next(iter(distinct[name]))}"
                parameter_box.addItem(f"{name}  [{tag}]", name)
            # Keeps the previous parameter, otherwise starts with the first one that varies
            index = parameter_box.findData(previous_parameter)
            if index < 0:
                index = next((i for i, n in enumerate(names) if len(distinct[n]) > 1), 0)
            parameter_box.setCurrentIndex(index)
            for name in available_plot_options(data['groups']):
                plot_box.addItem(name, name)  # the name is also the item's data, read by currentData()
            plot_box.setCurrentIndex(max(plot_box.findData(previous_plot), 0))

        for box in (parameter_box, plot_box):
            box.blockSignals(False)
        show_options(bool(indices))
        update_subspecies(keep_selection=True)
        if not indices:
            status.setText("No simulations selected.")
        elif not parameter_box.count() or not plot_box.count():
            status.setText("The selected simulations have no parameters in common or nothing to plot.")

    def update_subspecies(keep_selection=False):
        # Lists the subspecies of the chosen plot and shows only the menu rows it needs;
        # with keep_selection, the subspecies that were selected and are still listed stay selected
        plot = plot_box.currentData()
        previous = ({data['pairs'][subspecies_list.row(item)] for item in subspecies_list.selectedItems()}
                    if keep_selection else set())
        subspecies_list.blockSignals(True)  # the plot is updated once, at the end
        subspecies_list.clear()
        data['pairs'] = list_subspecies(plot_files(plot)[0], plot) if plot in SUBSPECIES_PLOTS else []
        subspecies_list.addItems([name for _, name in data['pairs']])
        for row, pair in enumerate(data['pairs']):
            if pair in previous:
                subspecies_list.item(row).setSelected(True)
        subspecies_list.blockSignals(False)
        menu.setRowVisible(subspecies_list, plot in SUBSPECIES_PLOTS)
        menu.setRowVisible(limit_box, plot == 'Particle Balance')
        update_time_rows()
        update_plot()

    def update_time_rows():
        # The animation options exist only for Final Densities; the time settings and the
        # save button only while animating
        densities = plot_box.currentData() == 'Final Densities'
        menu.setRowVisible(animate_check, densities)
        menu.setRowVisible(time_widget, densities and animate_check.isChecked())
        menu.setRowVisible(save_gif_button, densities and animate_check.isChecked())

    def stop_animation():
        # Stops the animation being played, if any (before the figure is cleared)
        anim = data['anim']
        if anim is not None and anim.event_source is not None:
            anim.pause()
        data['anim'] = None
        save_gif_button.setEnabled(False)

    def animate_densities(files, attrs, parameter, subspecies):
        # Builds the animation of the densities with the chosen time settings
        try:
            t_start, t_end = float(t_start_box.text()), float(t_end_box.text())
        except ValueError:
            return "Please give valid start and end times."
        instants, message = time_instants(t_start, t_end, n_frames_box.value())
        if instants is None:
            return message
        frames = read_time_densities(files, attrs, parameter, subspecies, instants)
        data['anim'], message = make_density_animation(figure, figure.add_subplot(), frames, instants,
                                                       parameter, subspecies)
        save_gif_button.setEnabled(data['anim'] is not None)
        return message

    def save_gif():
        # Asks where to save the animation (by default in the simulations folder) and saves it
        parameter = parameter_box.currentData()
        names = [name for _, name in (data['pairs'][subspecies_list.row(item)]
                                      for item in subspecies_list.selectedItems())]
        name = os.path.splitext(figure_filename('densities_anim', parameter, names))[0] + '.gif'
        default = os.path.join(data['path'], name)
        path, _ = QtWidgets.QFileDialog.getSaveFileName(window, "Save animation", default, "GIF (*.gif)")
        if not path or data['anim'] is None:  # cancelled
            return
        if not path.lower().endswith('.gif'):
            path += '.gif'
        status.setText("Saving the animation...")
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            save_animation(data['anim'], path)
            status.setText(f"Animation saved as {path}")
        except Exception as e:
            status.setText(f"Could not save the animation: {e}")
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()

    def update_plot():
        # Draws the chosen plot as soon as everything it needs is selected
        parameter, plot = parameter_box.currentData(), plot_box.currentData()
        subspecies = [data['pairs'][subspecies_list.row(item)] for item in subspecies_list.selectedItems()]
        message = note = ""
        stop_animation()
        figure.clear()
        if parameter is not None and plot is not None and (subspecies or plot == 'Reduced Field'):
            files, attrs, note = plot_files(plot)
            if plot == 'Final Densities' and animate_check.isChecked():
                message = animate_densities(files, attrs, parameter, subspecies)
            elif plot == 'Reduced Field':
                message = plot_reduced_field(figure.add_subplot(), files, attrs, parameter)
            elif plot in ('Final Densities', 'Populations'):
                message = plot_final_densities(figure.add_subplot(), files, attrs, parameter,
                                               subspecies, populations=(plot == 'Populations'))
            elif plot == 'Particle Balance':
                message = plot_particle_balance(figure.add_subplot(), files, attrs, parameter,
                                                subspecies, limit_box.value())
        status.setText(" ".join(text for text in (message, note) if text))
        figure.tight_layout()
        canvas.draw()
        toolbar.update()  # the "Home" button returns to this new plot

    def clear_plot():
        # Empties the figure and unselects the subspecies; changing any option draws again
        subspecies_list.blockSignals(True)  # so unselecting does not redraw
        subspecies_list.clearSelection()
        subspecies_list.blockSignals(False)
        stop_animation()
        figure.clear()
        status.setText("")
        canvas.draw()
        toolbar.update()

    folder_button.clicked.connect(open_folder)
    folder_tree.itemChanged.connect(lambda item, column: selection_timer.start())
    selection_timer.timeout.connect(apply_selection)
    all_button.clicked.connect(lambda: set_all_ticked(True))
    none_button.clicked.connect(lambda: set_all_ticked(False))
    clear_button.clicked.connect(clear_plot)
    save_gif_button.clicked.connect(save_gif)
    animate_check.toggled.connect(update_time_rows)
    animate_check.toggled.connect(update_plot)

    def time_edited(box):
        # Redraws when Enter is pressed or the box is left, only if its text was changed
        if box.isModified():
            box.setModified(False)
            update_plot()

    t_start_box.editingFinished.connect(lambda: time_edited(t_start_box))
    t_end_box.editingFinished.connect(lambda: time_edited(t_end_box))
    n_frames_box.valueChanged.connect(update_plot)
    parameter_box.currentIndexChanged.connect(update_plot)
    plot_box.currentIndexChanged.connect(lambda index: update_subspecies())
    subspecies_list.itemSelectionChanged.connect(update_plot)
    limit_box.valueChanged.connect(update_plot)

    show_options(False)  # only the first step until a folder is read
    menu.setRowVisible(tree_widget, False)
    update_subspecies()  # hides the subspecies and threshold rows until they are needed
    window.resize(1200, 700)
    window.show()
    sys.exit(app.exec())


def command_line():
    # Asks the same questions in the terminal (e.g. over ssh, where there is no display)
    # and then shows the figure in a window or saves it to a file
    # Wrapped to the terminal width (at most 90 columns; 80 if it cannot be detected, e.g. when piped)
    width = min(90, shutil.get_terminal_size().columns)
    print(textwrap.fill(INTRO_TEXT, width=width) + "\n")
    path_data = input("Path to simulations folder: ").strip().strip('"')
    if not os.path.isdir(path_data):
        print(f"Error: The specified path '{path_data}' does not exist.")
        sys.exit(1)

    readable_files, all_attrs, all_groups, _ = read_folder(path_data)
    if not readable_files:
        print(f"No chemistry simulations with .h5 files found in '{path_data}'.")
        sys.exit(1)
    print(f"Found {folder_summary(len(readable_files))}.")

    # --- Parameter for the x-axis ---
    names, distinct = common_attributes(all_attrs)
    if not names:
        print("There are no attributes common to all the files.")
        sys.exit(1)

    labels = []
    for name in names:
        n_values = len(distinct[name])
        tag = f"varies: {n_values} values" if n_values > 1 else f"constant = {next(iter(distinct[name]))}"
        labels.append(f"{name}  [{tag}]")

    parameter = names[ask_choice("Which parameter would you like to plot against (x axis)?", labels)]
    print(f"\nSelected parameter: {parameter}")

    # --- What to plot ---
    plot_options = available_plot_options(all_groups)
    if not plot_options:
        print("None of the plot options (Reduced Field, Final Densities, Populations, Particle Balance) "
              "is available in the files.")
        sys.exit(1)

    plot = plot_options[ask_choice("What do you want to plot (y axis)?", plot_options)]
    print(f"\nSelected plot: {plot}")

    # From here on, only the files with the data of the chosen plot are used
    readable_files, all_attrs, note = files_for_plot(plot, readable_files, all_attrs, all_groups)
    if note:
        print(note)

    selected_subspecies = []
    if plot in SUBSPECIES_PLOTS:
        pairs = list_subspecies(readable_files, plot)
        if not pairs:
            print(f"No subspecies found for {plot}.")
            sys.exit(1)

        selected_subspecies = ask_multiple_choice(
            "Which subspecies do you want to plot (in the same figure)?", pairs)
        print("\nSelected subspecies:")
        for parent, name in selected_subspecies:
            print(f"  {parent}: {name}")

    limit = None
    if plot == 'Particle Balance':
        limit = ask_limit()

    # --- Animation in time (Final Densities only) ---
    instants = None
    times = time_range(readable_files) if plot == 'Final Densities' else None
    if times is not None and ask_choice("Plot the final state or an animation in time?",
                                        ["Final state", "Animation in time"]) == 1:
        print(f"\nThe simulated times go from {times[0]:.3g} s to {times[1]:.3g} s.")
        while instants is None:
            t_start = ask_number("Start time (s): ", minimum=0)
            t_end = ask_number("End time (s): ", minimum=0)
            n_frames = ask_number("Number of frames: ", minimum=2, integer=True)
            instants, error = time_instants(t_start, t_end, n_frames)
            if instants is None:
                print(error)

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(10, 6) if plot == 'Particle Balance' else None)
    anim = None
    if instants is not None:
        frames = read_time_densities(readable_files, all_attrs, parameter, selected_subspecies, instants)
        anim, message = make_density_animation(fig, ax, frames, instants, parameter, selected_subspecies)
        if anim is None:
            print(message)
            sys.exit(1)
    elif plot == 'Reduced Field':
        message = plot_reduced_field(ax, readable_files, all_attrs, parameter)
    elif plot in ('Final Densities', 'Populations'):
        message = plot_final_densities(ax, readable_files, all_attrs, parameter, selected_subspecies,
                                       populations=(plot == 'Populations'))
    elif plot == 'Particle Balance':
        message = plot_particle_balance(ax, readable_files, all_attrs, parameter, selected_subspecies, limit)
    if message:
        print(message)
    if not ax.has_data():
        sys.exit(1)
    fig.tight_layout()

    # --- Show or save (without a display, e.g. over ssh, matplotlib uses 'agg' and can only save) ---
    if plt.get_backend().lower() == 'agg':
        print("\nNo display available: the figure will be saved.")
    elif ask_choice(f"Show the {'animation' if anim else 'figure'} or save it?", ["Show", "Save"]) == 0:
        plt.show()
        return

    kind = {'Reduced Field': 'field', 'Final Densities': 'densities', 'Populations': 'populations',
            'Particle Balance': 'balance'}[plot]
    if anim:
        kind += '_anim'
    extension = '.gif' if anim else '.png'
    default = os.path.splitext(figure_filename(kind, parameter, [name for _, name in selected_subspecies]))[0]
    default += extension
    name = input(f"\nFile name, in the simulations folder unless a full path is given "
                 f"(Enter for {default}): ").strip().strip('"') or default
    if not os.path.splitext(name)[1]:
        name += extension
    output_path = os.path.join(path_data, name)
    if anim:
        print("Saving the animation...")
        save_animation(anim, output_path)
    else:
        fig.savefig(output_path)
    print(f"\n{'Animation' if anim else 'Figure'} saved as: {output_path}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=INTRO_TEXT)
    parser.add_argument('--cli', action='store_true',
                        help="ask the questions in the terminal instead of opening a window (e.g. over ssh)")
    args = parser.parse_args()

    plt.rcParams.update({'savefig.dpi': 300, 'savefig.bbox': 'tight'})  # used by both save options
    if args.cli:
        command_line()
    else:
        gui()
