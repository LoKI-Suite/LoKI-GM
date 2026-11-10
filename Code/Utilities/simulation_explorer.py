"""

Interactive explorer for simulation results stored in HDF5 (.h5) files.

The script:
1) asks the user to select the folder where the simulations are stored
   (the .h5 files can be directly inside it or spread over subfolders);
2) walks through every folder and reads the attributes of each .h5 file;
3) lists the attributes common to all files and asks in function of which
   parameter the results should be represented;
4) asks what to plot (Reduced Field, Final Densities or Particle Balance) and,
   for densities/balance, which subspecies to plot in the same figure.

Run with:

python simulation_explorer.py

"""

__author__ = 'Pedro Pereira <pedro.t.pereira@tecnico.ulisboa.pt>'
__maintainer__ = 'Pedro Pereira <pedro.t.pereira@tecnico.ulisboa.pt>'
__date__ = '9 November 2026'
__version__ = 'Revision: 0.2.0'

import os
import re
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt


PLOT_OPTIONS = [
    ('Reduced Field', '/electronKinetics/reducedField'),
    ('Final Densities', '/chemistry/finalDensities'),
    ('Particle Balance', '/chemistry/particleBalance'),
]

FIELD_PATH = '/electronKinetics/reducedField'
DENSITIES_PATH = '/chemistry/finalDensities'
BALANCE_PATH = '/chemistry/particleBalance'
EXCLUDED_BALANCE_SUBSPECIES = {'N2(X)'}


def select_folder():
    # Opens a dialog to choose the simulations folder; falls back to typing the path
    path = ''
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        path = filedialog.askdirectory(title='Select the folder with the simulations')
        root.destroy()
    except Exception as e:
        print(f"Could not open the folder dialog ({e}).")

    if not path:
        path = input("Path to simulations folder: ").strip().strip('"')

    if not os.path.isdir(path):
        print(f"Error: The specified path '{path}' does not exist.")
        sys.exit(1)

    return path


def find_h5_files(path):
    # Walks through every folder inside path and collects the .h5 files
    h5_files = []
    for dirpath, dirnames, filenames in os.walk(path):
        for filename in filenames:
            if filename.endswith('.h5'):
                h5_files.append(os.path.join(dirpath, filename))

    if not h5_files:
        print(f"No .h5 files found in '{path}'.")
        sys.exit(1)

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


def plot_reduced_field(h5_files, all_attrs, parameter, path_fig):
    # Plots the reduced electric field E/N [Td] in function of the chosen parameter
    # (as in field_analyser of densities_field_cases.py)
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
        print("No reduced field information found in the dataset.")
        return

    order = np.argsort(x_values)
    x_values = np.asarray(x_values)[order]
    y_values = np.asarray(y_values)[order]

    fig, ax = plt.subplots()
    ax.plot(x_values, y_values, marker='o')
    ax.set_xlabel(parameter, fontsize=12)
    ax.set_ylabel("E/N (Td)", fontsize=12)
    ax.tick_params(labelsize=12)
    fig.tight_layout()

    show_figure(fig, path_fig, "field", parameter)


def figure_filename(kind, parameter, extra=()):
    # Builds e.g. "sim_field_Discharge_current.png": units between parenthesis are dropped
    # from the parameter, and any character not valid in a file name becomes "_"
    parameter_tag = re.sub(r'\([^)]*\)', '', parameter)
    parts = [re.sub(r'[^\w\-.]+', '_', part).strip('_')
             for part in ("sim", kind, parameter_tag, *extra)]
    return "_".join(part for part in parts if part) + ".png"


def show_figure(fig, path_fig, kind, parameter, extra=()):
    # Shows the figure in a Tk window (TkAgg, no extra install needed) with a side panel
    # to change the title and axis labels, the axis scales, show/hide the grid and save the figure
    import tkinter as tk
    from tkinter import ttk, filedialog
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    ax = fig.axes[0]

    root = tk.Tk()
    root.title("Simulation explorer")

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    panel = ttk.Frame(root, padding=10)
    panel.pack(side=tk.RIGHT, fill=tk.Y)
    status = tk.StringVar(master=root)

    # --- Labels (the figure updates while typing; an empty title means no title) ---
    labels_frame = ttk.LabelFrame(panel, text="Labels", padding=8)
    labels_frame.pack(fill=tk.X)
    label_setters = (('Title', ax.get_title(), ax.set_title),
                     ('x label', ax.get_xlabel(), ax.set_xlabel),
                     ('y label', ax.get_ylabel(), ax.set_ylabel))

    for row, (text, initial, setter) in enumerate(label_setters):
        ttk.Label(labels_frame, text=text).grid(row=2 * row, column=0, sticky='w')
        variable = tk.StringVar(master=root, value=initial)
        ttk.Entry(labels_frame, textvariable=variable, width=40).grid(row=2 * row + 1, column=0,
                                                                      sticky='we', pady=(0, 6))

        def apply_label(*_args, variable=variable, setter=setter):
            setter(variable.get(), fontsize=12)
            try:
                fig.tight_layout()
                canvas.draw()
                status.set("")
            except ValueError:  # e.g. unfinished LaTeX ($...$) while typing
                status.set("The label can't be rendered yet (incomplete LaTeX?).")

        variable.trace_add('write', apply_label)

    # --- Axis scales ---
    scales_frame = ttk.LabelFrame(panel, text="Axis scales", padding=8)
    scales_frame.pack(fill=tk.X, pady=(10, 0))
    x_scale = tk.StringVar(master=root, value=ax.get_xscale())
    y_scale = tk.StringVar(master=root, value=ax.get_yscale())

    def apply_scales(_event=None):
        ax.set_xscale(x_scale.get())
        ax.set_yscale(y_scale.get())
        ax.tick_params(which='both', labelsize=12)  # log scale adds minor tick labels
        # Values <= 0 cannot be drawn on a log axis (e.g. losses in the particle balance)
        hidden = []
        for axis, scale in (('x', x_scale.get()), ('y', y_scale.get())):
            values = [line.get_xdata() if axis == 'x' else line.get_ydata() for line in ax.get_lines()]
            if scale == 'log' and any(np.any(np.asarray(v) <= 0) for v in values):
                hidden.append(axis)
        status.set(f"Values <= 0 on the {' and '.join(hidden)} axis are not shown in log scale."
                   if hidden else "")
        fig.tight_layout()
        canvas.draw_idle()

    for row, (axis, variable) in enumerate((('x axis', x_scale), ('y axis', y_scale))):
        ttk.Label(scales_frame, text=axis).grid(row=row, column=0, sticky='w', padx=(0, 8), pady=2)
        box = ttk.Combobox(scales_frame, textvariable=variable, values=('linear', 'log'),
                           state='readonly', width=8)
        box.grid(row=row, column=1, sticky='w', pady=2)
        box.bind('<<ComboboxSelected>>', apply_scales)

    # --- Grid ---
    grid_frame = ttk.LabelFrame(panel, text="Grid", padding=8)
    grid_frame.pack(fill=tk.X, pady=(10, 0))
    show_grid = tk.BooleanVar(master=root, value=False)

    def apply_grid():
        ax.grid(show_grid.get())
        canvas.draw_idle()

    ttk.Checkbutton(grid_frame, text="Show grid", variable=show_grid,
                    command=apply_grid).pack(anchor='w')

    # --- Save ---
    save_frame = ttk.LabelFrame(panel, text="Save", padding=8)
    save_frame.pack(fill=tk.X, pady=(10, 0))
    file_name = tk.StringVar(master=root, value=figure_filename(kind, parameter, extra))
    dpi = tk.StringVar(master=root, value='300')
    last_folder = [path_fig]

    ttk.Label(save_frame, text="File name").grid(row=0, column=0, columnspan=2, sticky='w')
    ttk.Entry(save_frame, textvariable=file_name, width=40).grid(row=1, column=0, columnspan=2,
                                                                  sticky='we', pady=(0, 6))
    ttk.Label(save_frame, text="DPI").grid(row=2, column=0, sticky='w', padx=(0, 8))
    ttk.Entry(save_frame, textvariable=dpi, width=8).grid(row=2, column=1, sticky='w')

    def save():
        name = file_name.get().strip()
        if not name:
            status.set("Write a file name.")
            return
        try:
            dpi_value = float(dpi.get())
            if dpi_value <= 0:
                raise ValueError
        except ValueError:
            status.set("DPI must be a positive number.")
            return
        if not os.path.splitext(name)[1]:
            name += '.png'
        # Window to choose where to save; it starts in the last folder used
        # (the simulations folder the first time) with the name written above
        output_path = filedialog.asksaveasfilename(
            parent=root, title="Save figure", initialdir=last_folder[0], initialfile=name,
            defaultextension='.png',
            filetypes=[('PNG', '*.png'), ('PDF', '*.pdf'), ('SVG', '*.svg'), ('All files', '*.*')])
        if not output_path:  # cancelled
            return
        try:
            fig.savefig(output_path, dpi=dpi_value, bbox_inches='tight')
        except Exception as e:
            status.set(f"Could not save the figure ({e}).")
            return
        last_folder[0] = os.path.dirname(output_path)
        file_name.set(os.path.basename(output_path))
        status.set(f"Saved: {output_path}")
        print(f"\nFigure saved as: {output_path}")

    ttk.Button(save_frame, text="Save", command=save).grid(row=3, column=0, columnspan=2,
                                                           sticky='we', pady=(8, 0))

    ttk.Label(panel, textvariable=status, wraplength=260).pack(fill=tk.X, pady=(10, 0))

    def on_close():
        root.quit()
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", on_close)
    root.mainloop()
    plt.close(fig)


def read_density(f, parent, name):
    # Absolute density of a subspecies, found by name (None if not in this file)
    if parent not in f[DENSITIES_PATH]:
        return None
    dataset = f[DENSITIES_PATH][parent]
    if dataset.ndim < 2 or dataset.shape[-1] == 0:
        return None
    for i in range(dataset.shape[1]):
        row = dataset[0, i]
        if row['Species'].decode('utf-8') == name:
            return float(row['Abs.Density'])
    return None


def plot_final_densities(h5_files, all_attrs, parameter, subspecies, path_fig):
    # Plots the absolute density of several subspecies in function of the chosen parameter,
    # all in the same figure (as densities_analyser with -hold in densities_field_cases.py)
    data = {pair: ([], []) for pair in subspecies}
    for file_path, attrs in zip(h5_files, all_attrs):
        with h5py.File(file_path, 'r') as f:
            for parent, name in subspecies:
                value = read_density(f, parent, name)
                if value is None:
                    continue
                data[(parent, name)][0].append(attrs[parameter])
                data[(parent, name)][1].append(value)

    if not any(x for x, _ in data.values()):
        print("No density information found in the dataset for the selected subspecies.")
        return

    fig, ax = plt.subplots()
    for (parent, name), (x_values, y_values) in data.items():
        if not x_values:
            print(f"No density information found for {name}.")
            continue
        order = np.argsort(x_values)
        ax.plot(np.asarray(x_values)[order], np.asarray(y_values)[order], marker='o', label=name)

    # Wall species are surface densities (m-2), the others volume densities (m-3)
    units = {'m$^{-2}$' if parent.lower().startswith('wall') else 'm$^{-3}$'
             for parent, _ in subspecies}
    ax.set_xlabel(parameter, fontsize=12)
    ax.set_ylabel(f"Absolute density ({' / '.join(sorted(units))})", fontsize=12)
    ax.tick_params(labelsize=12)
    ax.yaxis.get_offset_text().set_fontsize(12)  # the "1e22" multiplier
    ax.legend(fontsize=12)
    fig.tight_layout()

    show_figure(fig, path_fig, "densities", parameter, [name for _, name in subspecies])


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


def plot_particle_balance(h5_files, all_attrs, parameter, subspecies, limit, path_fig):
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
        print(f"No reactions above {limit}% found for the selected subspecies.")
        return

    # The default palette has 10 colours; with more lines use tab20 so they don't repeat
    n_default = len(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    colors = plt.cm.tab20.colors if len(data) > n_default else None

    fig, ax = plt.subplots(figsize=(10, 6))
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
    ax.legend()  # placed automatically inside the plot, default font size
    fig.tight_layout()

    show_figure(fig, path_fig, "balance", parameter, [name for _, name in subspecies])


if __name__ == "__main__":

    path_data = select_folder()
    print(f"Simulations folder: {path_data}")

    h5_files = find_h5_files(path_data)
    print(f"Found {len(h5_files)} .h5 file(s).")

    readable_files = []
    all_attrs = []
    all_groups = []
    for file_path in h5_files:
        try:
            attrs, groups = read_file_info(file_path)
        except OSError as e:
            print(f"Warning: could not read '{file_path}' ({e}). Skipping.")
            continue
        readable_files.append(file_path)
        all_attrs.append(attrs)
        all_groups.append(groups)

    if not all_attrs:
        print("No readable .h5 files found.")
        sys.exit(1)

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
        print("None of the plot options (Reduced Field, Final Densities, Particle Balance) "
              "is available in all the files.")
        sys.exit(1)

    plot = plot_options[ask_choice("What do you want to plot?", plot_options)]
    print(f"\nSelected plot: {plot}")

    selected_subspecies = []
    if plot in ('Final Densities', 'Particle Balance'):
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
    if plot == 'Reduced Field':
        plot_reduced_field(readable_files, all_attrs, parameter, path_data)
    elif plot == 'Final Densities':
        plot_final_densities(readable_files, all_attrs, parameter, selected_subspecies, path_data)
    elif plot == 'Particle Balance':
        plot_particle_balance(readable_files, all_attrs, parameter, selected_subspecies, limit, path_data)
