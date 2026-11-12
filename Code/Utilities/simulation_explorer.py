"""

Interactive explorer for simulation results stored in HDF5 (.h5) files.

The window has a menu on the left and the figure on the right. In the menu:
1) "Open folder..." selects the folder where the simulations are stored
   (the .h5 files can be directly inside it or spread over subfolders);
2) "Parameter" chooses, among the attributes common to all files, in function
   of which one the results are represented;
3) "Plot" chooses what to plot (Reduced Field, Final Densities, Populations or
   Particle Balance) and, for densities/populations/balance, which subspecies to plot in the same figure.
The figure is drawn as soon as everything is selected; its toolbar allows to zoom,
edit the title, labels and axis scales, and save it.

The same questions can also be answered in the terminal (e.g. over ssh, without a
display); the figure is then shown in a window or saved to a file.

Run with:

python simulation_explorer.py          (window)
python simulation_explorer.py --cli    (terminal)

"""

__author__ = 'Pedro Pereira <pedro.t.pereira@tecnico.ulisboa.pt>'
__maintainer__ = 'Pedro Pereira <pedro.t.pereira@tecnico.ulisboa.pt>'
__date__ = '9 November 2026'
__version__ = 'Revision: 0.2.0'

import argparse
import os
import re
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt


PLOT_OPTIONS = [
    ('Reduced Field', '/electronKinetics/reducedField'),
    ('Final Densities', '/chemistry/finalDensities'),
    ('Populations', '/chemistry/finalDensities'),
    ('Particle Balance', '/chemistry/particleBalance'),
]

FIELD_PATH = '/electronKinetics/reducedField'
DENSITIES_PATH = '/chemistry/finalDensities'
BALANCE_PATH = '/chemistry/particleBalance'
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
    # Reads the root attributes of an .h5 file and which plot groups it contains
    with h5py.File(file_path, 'r') as f:
        attrs = {name: to_python_value(f.attrs[name]) for name in f.attrs}
        groups = {name for name, path in PLOT_OPTIONS if path in f}
    return attrs, groups


def read_folder(path):
    # Reads the info of every .h5 file inside path, skipping the ones that cannot be opened
    files, all_attrs, all_groups = [], [], []
    for file_path in find_h5_files(path):
        try:
            attrs, groups = read_file_info(file_path)
        except OSError as e:
            print(f"Warning: could not read '{file_path}' ({e}). Skipping.")
            continue
        files.append(file_path)
        all_attrs.append(attrs)
        all_groups.append(groups)
    return files, all_attrs, all_groups


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
    # Keeps only the plot options whose group exists in every file
    common = set.intersection(*all_groups)
    return [name for name, _ in PLOT_OPTIONS if name in common]


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
    if populations:
        ax.set_ylabel("Population", fontsize=12)
    else:
        # Wall species are surface densities (m-2), the others volume densities (m-3)
        units = {'m$^{-2}$' if parent.lower().startswith('wall') else 'm$^{-3}$'
                 for parent, _ in subspecies}
        ax.set_ylabel(f"Absolute density ({' / '.join(sorted(units))})", fontsize=12)
    ax.tick_params(labelsize=12)
    ax.yaxis.get_offset_text().set_fontsize(12)  # the "1e22" multiplier
    ax.legend(fontsize=12)
    return f"No {quantity} information found for {', '.join(missing)}." if missing else ""


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
    from matplotlib.backends.qt_compat import QtWidgets
    from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
    from matplotlib.figure import Figure

    app = QtWidgets.QApplication(sys.argv)

    window = QtWidgets.QWidget(windowTitle="Simulation explorer")
    folder_button = QtWidgets.QPushButton("Open folder...")
    folder_label = QtWidgets.QLabel("No folder selected")
    parameter_box = QtWidgets.QComboBox()  # long names are cut with "..." but shown complete when opened
    parameter_box.setSizeAdjustPolicy(QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
    plot_box = QtWidgets.QComboBox()
    subspecies_list = QtWidgets.QListWidget()  # several can be selected, one click each
    subspecies_list.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)
    limit_box = QtWidgets.QDoubleSpinBox(maximum=100, value=5, suffix=' %')
    clear_button = QtWidgets.QPushButton("Clear plot")
    status = QtWidgets.QLabel(wordWrap=True)

    menu_widget = QtWidgets.QWidget()
    menu_widget.setFixedWidth(300)  # so it does not change when rows are shown or hidden
    menu = QtWidgets.QFormLayout(menu_widget)
    menu.setContentsMargins(0, 0, 0, 0)
    menu.setRowWrapPolicy(QtWidgets.QFormLayout.RowWrapPolicy.WrapAllRows)  # labels above the widgets
    menu.addRow(folder_button)
    menu.addRow(folder_label)
    menu.addRow("Parameter (x axis)", parameter_box)
    menu.addRow("Plot", plot_box)
    menu.addRow("Subspecies", subspecies_list)
    menu.addRow("Show reactions above", limit_box)
    menu.addRow(clear_button)
    menu.addRow(status)

    figure = Figure()
    canvas = FigureCanvasQTAgg(figure)
    toolbar = NavigationToolbar2QT(canvas, window)
    plot_area = QtWidgets.QVBoxLayout()
    plot_area.addWidget(toolbar)
    plot_area.addWidget(canvas)

    layout = QtWidgets.QHBoxLayout(window)
    layout.addWidget(menu_widget)
    layout.addLayout(plot_area, stretch=1)

    data = {'files': [], 'attrs': [], 'pairs': []}

    def open_folder():
        # Reads the chosen folder and fills the menus with its parameters and plot options
        path = QtWidgets.QFileDialog.getExistingDirectory(window, "Select the folder with the simulations")
        if not path:  # cancelled
            return
        files, all_attrs, all_groups = read_folder(path)
        if not files:
            status.setText(f"No readable .h5 files found in '{path}'.")
            return
        names, distinct = common_attributes(all_attrs)
        data.update(files=files, attrs=all_attrs)
        folder_label.setText(f"{os.path.basename(path)}: {len(files)} .h5 file(s)")
        plt.rcParams['savefig.directory'] = path  # where the save dialog starts

        plot_box.clear()  # first, so nothing is plotted while the parameters are refilled
        parameter_box.clear()
        for name in names:
            n_values = len(distinct[name])
            tag = f"varies: {n_values} values" if n_values > 1 else f"constant = {next(iter(distinct[name]))}"
            parameter_box.addItem(f"{name}  [{tag}]", name)
        # Starts with the first parameter that varies
        parameter_box.setCurrentIndex(next((i for i, n in enumerate(names) if len(distinct[n]) > 1), 0))
        plot_box.addItems(available_plot_options(all_groups))
        if not names or not plot_box.count():
            status.setText("The files have no parameters or plot options in common.")

    def update_subspecies():
        # Lists the subspecies of the chosen plot and shows only the menu rows it needs
        plot = plot_box.currentText()
        subspecies_list.clear()
        data['pairs'] = list_subspecies(data['files'], plot) if plot in SUBSPECIES_PLOTS else []
        subspecies_list.addItems([name for _, name in data['pairs']])
        menu.setRowVisible(subspecies_list, plot in SUBSPECIES_PLOTS)
        menu.setRowVisible(limit_box, plot == 'Particle Balance')
        update_plot()

    def update_plot():
        # Draws the chosen plot as soon as everything it needs is selected
        parameter, plot = parameter_box.currentData(), plot_box.currentText()
        subspecies = [data['pairs'][subspecies_list.row(item)] for item in subspecies_list.selectedItems()]
        message = ""
        figure.clear()
        if parameter is not None and plot == 'Reduced Field':
            message = plot_reduced_field(figure.add_subplot(), data['files'], data['attrs'], parameter)
        elif parameter is not None and subspecies and plot in ('Final Densities', 'Populations'):
            message = plot_final_densities(figure.add_subplot(), data['files'], data['attrs'], parameter,
                                           subspecies, populations=(plot == 'Populations'))
        elif parameter is not None and subspecies and plot == 'Particle Balance':
            message = plot_particle_balance(figure.add_subplot(), data['files'], data['attrs'], parameter,
                                            subspecies, limit_box.value())
        status.setText(message)
        figure.tight_layout()
        canvas.draw()
        toolbar.update()  # the "Home" button returns to this new plot

    def clear_plot():
        # Empties the figure and unselects the subspecies; changing any option draws again
        subspecies_list.blockSignals(True)  # so unselecting does not redraw
        subspecies_list.clearSelection()
        subspecies_list.blockSignals(False)
        figure.clear()
        status.setText("")
        canvas.draw()
        toolbar.update()

    folder_button.clicked.connect(open_folder)
    clear_button.clicked.connect(clear_plot)
    parameter_box.currentIndexChanged.connect(update_plot)
    plot_box.currentIndexChanged.connect(update_subspecies)
    subspecies_list.itemSelectionChanged.connect(update_plot)
    limit_box.valueChanged.connect(update_plot)

    update_subspecies()  # hides the subspecies and threshold rows until they are needed
    window.resize(1200, 700)
    window.show()
    sys.exit(app.exec())


def command_line():
    # Asks the same questions in the terminal (e.g. over ssh, where there is no display)
    # and then shows the figure in a window or saves it to a file
    path_data = input("Path to simulations folder: ").strip().strip('"')
    if not os.path.isdir(path_data):
        print(f"Error: The specified path '{path_data}' does not exist.")
        sys.exit(1)

    readable_files, all_attrs, all_groups = read_folder(path_data)
    if not readable_files:
        print("No readable .h5 files found.")
        sys.exit(1)
    print(f"Found {len(readable_files)} .h5 file(s).")

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

    parameter = names[ask_choice("Should we represent in function of which parameter?", labels)]
    print(f"\nSelected parameter: {parameter}")

    # --- What to plot ---
    plot_options = available_plot_options(all_groups)
    if not plot_options:
        print("None of the plot options (Reduced Field, Final Densities, Populations, Particle Balance) "
              "is available in all the files.")
        sys.exit(1)

    plot = plot_options[ask_choice("What do you want to plot?", plot_options)]
    print(f"\nSelected plot: {plot}")

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

    # --- Plotting ---
    fig, ax = plt.subplots(figsize=(10, 6) if plot == 'Particle Balance' else None)
    if plot == 'Reduced Field':
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
    elif ask_choice("Show the figure or save it?", ["Show", "Save"]) == 0:
        plt.show()
        return

    kind = {'Reduced Field': 'field', 'Final Densities': 'densities', 'Populations': 'populations',
            'Particle Balance': 'balance'}[plot]
    default = figure_filename(kind, parameter, [name for _, name in selected_subspecies])
    name = input(f"\nFile name, in the simulations folder unless a full path is given "
                 f"(Enter for {default}): ").strip().strip('"') or default
    if not os.path.splitext(name)[1]:
        name += '.png'
    output_path = os.path.join(path_data, name)
    fig.savefig(output_path)
    print(f"\nFigure saved as: {output_path}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Explorer for simulation results stored in HDF5 (.h5) files.")
    parser.add_argument('--cli', action='store_true',
                        help="ask the questions in the terminal instead of opening a window (e.g. over ssh)")
    args = parser.parse_args()

    plt.rcParams.update({'savefig.dpi': 300, 'savefig.bbox': 'tight'})  # used by both save options
    if args.cli:
        command_line()
    else:
        gui()
