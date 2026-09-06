# -*- coding: utf-8 -*-
"""
AllyMind 10.2 FINAL – контекстное меню, Control Zero fix, ядерные изомеры (встроенные),
двойные мутации (double mutant), статистика анализов, таймеры, поддержка PDB _A.
Confidential. Valentin Lebedkin, 11.07.2026.
"""
import sys, os, json, subprocess, webbrowser, tempfile, tkinter as tk
from tkinter import ttk, messagebox, filedialog, font as tkfont
import sqlite3, io, datetime, time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from mpl_toolkits.mplot3d import Axes3D
import allymind_core
from allymind_core import AA_PROPS, determine_strategy
import cascade
import pdb_parser

try:
    from Bio.PDB import PDBParser
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

def resource_path(relative_path):
    try: return os.path.join(sys._MEIPASS, relative_path)
    except: return os.path.abspath(relative_path)

APP_VERSION = "10.2 FINAL – stats + double mutant + timer + PDB _A"
AUTHOR = "Valentin Lebedkin"
COPYRIGHT = f"© {AUTHOR}, 2026. All rights reserved.\nConfidential."

def check_trap_old(sequence, pos, ss_map=None):
    if ss_map and 0 <= pos < len(ss_map):
        J_local, _ = cascade.secondary_structure_J(sequence, pos, ss_map)
    else:
        J_local, _ = cascade.secondary_structure_J(sequence, pos, None)
    dom_rad = 10
    N = len(sequence)
    start = max(0, pos - dom_rad)
    end = min(N, pos + dom_rad + 1)
    dom_len = end - start
    healthy = N - dom_len
    defect = {'nh': healthy, 'nd': dom_len, 'J': J_local, 'N': N}
    hole   = {'nh': healthy, 'nd': 0, 'J': cascade.J_NAT, 'N': N}
    F_def = cascade.dimer_free_energy(defect, defect, cascade.J_NAT, dom_len)
    F_mix = cascade.dimer_free_energy(defect, hole, cascade.J_NAT, dom_len)
    return F_def < F_mix, F_mix - F_def

def check_trap_new(sequence, pos, ref, alt, ss_map=None):
    mut_code = f"{ref}{pos+1}{alt}"
    report = allymind_core.analyze_mutation(sequence, mut_code, ss_map=ss_map)
    return report['trap'], report['depth'], cascade.is_trap_exposed(sequence, pos, ref, alt, ss_map=ss_map)

class AllyMindApp:
    def __init__(self, root):
        self.root = root
        self.root.title(f"AllyMind {APP_VERSION}")
        self.root.geometry("1400x900")
        self.root.configure(bg="#2d2d2d")
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.root.option_add("*TEntry*insertColor", "#ffffff")
        self.set_english_keyboard_layout()
        if not self.show_warning(): sys.exit(0)
        self.font_default = tkfont.Font(family="Segoe UI", size=10)
        self.font_bold   = tkfont.Font(family="Segoe UI", size=10, weight="bold")
        self.font_title  = tkfont.Font(family="Segoe UI", size=18, weight="bold")
        self.style = ttk.Style()
        self.style.theme_use("clam")
        self.style.configure(".", background="#2d2d2d", foreground="#eff0f1", font=self.font_default)
        self.style.configure("TNotebook", background="#2d2d2d", borderwidth=0)
        self.style.configure("TNotebook.Tab", background="#3e3e3e", foreground="#eff0f1", padding=[12,6], font=self.font_bold)
        self.style.map("TNotebook.Tab", background=[("selected","#2d2d2d")], foreground=[("selected","#3daee9")])
        self.style.configure("TButton", background="#3daee9", foreground="#ffffff", borderwidth=0, relief="flat", padding=[10,5])
        self.style.map("TButton", background=[("active","#2980b9")], foreground=[("active","#ffffff")])
        self.style.configure("TLabel", background="#2d2d2d", foreground="#eff0f1")
        self.style.configure("TEntry", fieldbackground="#3e3e3e", foreground="#eff0f1", borderwidth=1, relief="solid")
        self.style.configure("Treeview", background="#3e3e3e", foreground="#eff0f1", fieldbackground="#3e3e3e")
        self.style.configure("Treeview.Heading", background="#2d2d2d", foreground="#3daee9", font=self.font_bold)
        self.style.configure("TFrame", background="#2d2d2d")
        icon_path = resource_path("allymind_logo.ico")
        if os.path.exists(icon_path): self.root.iconbitmap(icon_path)
        self.sequence = ""
        self.gene_symbol = ""
        self.uniprot_id = ""
        self.current_mutation = ""
        self.ss_map = None
        self.pdb_path = None
        self.metal_sites = []
        self.metal_coordinated_residues = set()
        self.dimer_interface_residues = set()
        self.analysis_count = {}
        self._last_fig = None
        self._last_scatter_data = None
        self._last_full_reports = []
        self._pan_active = False
        self._rot_active = False
        self._pan_start_x = 0
        self._pan_start_y = 0
        self.build_ui()
        self.entry_id.focus_force()
        self.root.update_idletasks()

    def set_english_keyboard_layout(self):
        try:
            import ctypes
            ctypes.windll.user32.LoadKeyboardLayoutW("00000409", 0x00000001)
        except: pass

    def show_warning(self):
        return messagebox.askokcancel(
            "CONFIDENTIAL — ALL RIGHTS RESERVED",
            "AllyMind 10.2 FINAL – stats + double mutant + timer + PDB _A\n\n"
            "© Valentin Lebedkin, 2026. All rights reserved. Confidential.\n"
            "© Валентин Лебедкин, 2026. Все права защищены. Конфиденциально.\n\n"
            "This software contains trade secrets.\n"
            "Это программное обеспечение содержит коммерческую тайну.\n\n"
            "Do not distribute, decompile, or modify.\n"
            "Не распространять, не декомпилировать и не изменять.\n\n"
            "Continue? / Продолжить?",
            icon="warning"
        )

    def on_close(self):
        os.system("taskkill /F /IM python.exe /T >nul 2>&1")
        self.root.destroy()
        os._exit(0)

    def _increment_counter(self, gene=None):
        g = gene or self.gene_symbol
        if g:
            self.analysis_count[g] = self.analysis_count.get(g, 0) + 1

    def show_analysis_stats(self):
        if not self.analysis_count:
            messagebox.showinfo("Stats", "No analyses recorded yet.")
            return
        lines = [f"{gene}: {cnt} analysis(es)" for gene, cnt in self.analysis_count.items()]
        messagebox.showinfo("Analysis Statistics", "\n".join(lines))

    def analyze_double(self):
        start_time = time.time()
        mut1 = self.entry_manual.get().strip()
        mut2 = self.entry_manual2.get().strip()
        if not mut2 and ' ' in mut1:
            parts = mut1.split()
            if len(parts) >= 2:
                mut1, mut2 = parts[0], parts[1]
        if not mut1 or not mut2 or len(mut1) < 3 or len(mut2) < 3:
            messagebox.showwarning("Double", "Enter both mutations (e.g., A4V I113T)."); return
        if not self.sequence:
            messagebox.showwarning("Double", "Load a protein first."); return

        wt1, pos1, alt1 = mut1[0].upper(), int(mut1[1:-1]) - 1, mut1[-1].upper()
        wt2, pos2, alt2 = mut2[0].upper(), int(mut2[1:-1]) - 1, mut2[-1].upper()
        r1 = allymind_core.analyze_mutation(self.sequence, mut1.upper(), ss_map=self.ss_map)
        r2 = allymind_core.analyze_mutation(self.sequence, mut2.upper(), ss_map=self.ss_map)
        trap_old_1, depth_old_1 = check_trap_old(self.sequence, pos1, self.ss_map)
        trap_old_2, depth_old_2 = check_trap_old(self.sequence, pos2, self.ss_map)
        depth_new_1 = r1['depth'] if wt1 != alt1 else depth_old_1
        depth_new_2 = r2['depth'] if wt2 != alt2 else depth_old_2
        exposed_1 = cascade.is_trap_exposed(self.sequence, pos1, wt1, alt1, self.ss_map)
        exposed_2 = cascade.is_trap_exposed(self.sequence, pos2, wt2, alt2, self.ss_map)
        combined_depth = depth_new_1 + depth_new_2
        combined_ddG = (r1.get('kinetics', {}).get('delta_G', 0) or 0) + (r2.get('kinetics', {}).get('delta_G', 0) or 0)
        combined_trap = r1['trap'] or r2['trap']
        combined_exposed = exposed_1 or exposed_2
        classification = "open" if combined_exposed else "hidden"
        thr = max(r1.get('therapeutic_threshold', 0) or 0, r2.get('therapeutic_threshold', 0) or 0)
        dG = combined_ddG * 0.5
        ss_1 = r1.get('secondary_structure', 'loop') or 'loop'
        ss_2 = r2.get('secondary_structure', 'loop') or 'loop'
        comp_1 = r1.get('chemical_components', {})
        comp_2 = r2.get('chemical_components', {})

        strat1 = determine_strategy(
            mut1.upper(), depth_old_1, depth_new_1, "open" if exposed_1 else "hidden", thr,
            self.uniprot_id, self.entry_pdb.get().strip().upper(),
            self.sequence, self.ss_map, self.pdb_path,
            self.metal_coordinated_residues,
            dimer_interface_residues=self.dimer_interface_residues)
        strat2 = determine_strategy(
            mut2.upper(), depth_old_2, depth_new_2, "open" if exposed_2 else "hidden", thr,
            self.uniprot_id, self.entry_pdb.get().strip().upper(),
            self.sequence, self.ss_map, self.pdb_path,
            self.metal_coordinated_residues,
            dimer_interface_residues=self.dimer_interface_residues)

        if "Small molecule" in strat1 or "Small molecule" in strat2:
            combined_strategy = "Small molecule (combined)"
        elif "Chaperone" in strat1 or "Chaperone" in strat2:
            combined_strategy = "Chaperone (combined)"
        elif "Monitor" in strat1 and "Monitor" in strat2:
            combined_strategy = "Monitor"
        else:
            combined_strategy = f"{strat1} + {strat2}"

        self.tree.delete(*self.tree.get_children())
        self.tree.insert("", "end", values=(pos1+1, f"{wt1}\u2192{alt1}", f"{depth_old_1:.3f}", f"{depth_new_1:.3f}", ss_1, "open" if exposed_1 else "hidden", f"{r1.get('therapeutic_threshold', 0) or 0:.4f}", f"{(r1.get('kinetics', {}) or {}).get('delta_G', 0) or 0:.4f}", f"{r1.get('delta_delta_G', 0) or 0:.4f}", f"{comp_1.get('delta_hyd', 0):.4f}", f"{comp_1.get('delta_vol', 0):.4f}", f"{comp_1.get('ss_bonus', 0):.4f}", f"{comp_1.get('helix_penalty', 0):.4f}", strat1))
        self.tree.insert("", "end", values=(pos2+1, f"{wt2}\u2192{alt2}", f"{depth_old_2:.3f}", f"{depth_new_2:.3f}", ss_2, "open" if exposed_2 else "hidden", f"{r2.get('therapeutic_threshold', 0) or 0:.4f}", f"{(r2.get('kinetics', {}) or {}).get('delta_G', 0) or 0:.4f}", f"{r2.get('delta_delta_G', 0) or 0:.4f}", f"{comp_2.get('delta_hyd', 0):.4f}", f"{comp_2.get('delta_vol', 0):.4f}", f"{comp_2.get('ss_bonus', 0):.4f}", f"{comp_2.get('helix_penalty', 0):.4f}", strat2))
        self.tree.insert("", "end", values=("\u2014", f"{mut1}+{mut2}", f"{depth_old_1 + depth_old_2:.3f}", f"{combined_depth:.3f}", f"{ss_1}+{ss_2}", classification, f"{thr:.4f}", f"{dG:.4f}", f"{combined_ddG:.4f}", "\u2014", "\u2014", "\u2014", "\u2014", combined_strategy))

        self.current_mutation = f"{mut1}+{mut2}"
        self.plot_landscape_for_mutation(pos1)
        self._update_3d_structure()
        self._plot_stpc_graph()
        self._increment_counter()
        elapsed = (time.time() - start_time) * 1000
        self.set_status(f"Double: {mut1}+{mut2} | Trap: {'yes' if combined_trap else 'no'} | {elapsed:.1f} ms")

    def build_ui(self):
        notebook = ttk.Notebook(self.root)
        notebook.pack(fill="both", expand=True, padx=10, pady=10)
        self.tab_analysis = ttk.Frame(notebook); notebook.add(self.tab_analysis, text="Analysis"); self.build_analysis_tab()
        self.tab_graph    = ttk.Frame(notebook); notebook.add(self.tab_graph, text="Graph"); self.build_graph_tab()
        self.tab_3d       = ttk.Frame(notebook); notebook.add(self.tab_3d, text="3D Structure"); self.build_3d_tab()
        self.tab_stpc     = ttk.Frame(notebook); notebook.add(self.tab_stpc, text="STPC Model"); self.build_stpc_tab()
        self.tab_db       = ttk.Frame(notebook); notebook.add(self.tab_db, text="Database"); self.build_db_tab()
        self.tab_reports  = ttk.Frame(notebook); notebook.add(self.tab_reports, text="Reports"); self.build_reports_tab()
        self.tab_about    = ttk.Frame(notebook); notebook.add(self.tab_about, text="About"); self.build_about_tab()
        self.tab_isomer   = ttk.Frame(notebook); notebook.add(self.tab_isomer, text="Ядерные изомеры"); self.build_isomer_tab()
        status_frame = tk.Frame(self.root, bg="#3e3e3e"); status_frame.pack(fill="x", side="bottom")
        self.status = tk.Label(status_frame, text="Ready", anchor="w", bg="#3e3e3e", fg="#eff0f1", font=self.font_default)
        self.status.pack(side="left", padx=10, pady=5)
        tk.Label(status_frame, text=f"© {AUTHOR}, 2026", anchor="e", bg="#3e3e3e", fg="#888888").pack(side="right", padx=10, pady=5)
        self.notebook = notebook
        self.notebook.bind('<<NotebookTabChanged>>', self.on_tab_changed)

    def _add_context_menu(self, widget):
        menu = tk.Menu(widget, tearoff=0, bg="#3e3e3e", fg="#eff0f1")
        menu.add_command(label="Cut", command=lambda: widget.event_generate('<<Cut>>'))
        menu.add_command(label="Copy", command=lambda: widget.event_generate('<<Copy>>'))
        menu.add_command(label="Paste", command=lambda: widget.event_generate('<<Paste>>'))
        menu.add_command(label="Select All", command=lambda: widget.event_generate('<<SelectAll>>'))
        def show_menu(event):
            menu.tk_popup(event.x_root, event.y_root)
        widget.bind("<Button-3>", show_menu)
        widget.bind("<Control-v>", lambda e: widget.event_generate('<<Paste>>'))
        widget.bind("<Control-V>", lambda e: widget.event_generate('<<Paste>>'))
        widget.bind("<Control-c>", lambda e: widget.event_generate('<<Copy>>'))
        widget.bind("<Control-C>", lambda e: widget.event_generate('<<Copy>>'))

    def build_analysis_tab(self):
        frame = ttk.Frame(self.tab_analysis, padding=20); frame.pack(fill="both", expand=True)
        tk.Label(frame, text="AllyMind", font=self.font_title, fg="#3daee9", bg="#2d2d2d").grid(row=0, column=0, columnspan=14, pady=(0,15), sticky="w")
        ttk.Label(frame, text="UniProt ID or gene symbol:").grid(row=1, column=0, sticky="w", pady=5)
        self.entry_id = tk.Entry(frame, width=25, bg='#3e3e3e', fg='#eff0f1', insertbackground='white', relief='solid', bd=1)
        self.entry_id.grid(row=1, column=1, sticky="w", padx=10)
        self._add_context_menu(self.entry_id)
        ttk.Button(frame, text="Load", command=self.load_protein).grid(row=1, column=2, padx=5)
        ttk.Button(frame, text="FASTA...", command=self.load_fasta).grid(row=1, column=3, padx=5)
        ttk.Label(frame, text="PDB ID (optional):").grid(row=2, column=0, sticky="w", pady=5)
        self.entry_pdb = tk.Entry(frame, width=10, bg='#3e3e3e', fg='#eff0f1', insertbackground='white', relief='solid', bd=1)
        self.entry_pdb.grid(row=2, column=1, sticky="w", padx=10)
        self._add_context_menu(self.entry_pdb)
        ttk.Button(frame, text="Load PDB", command=self.load_pdb).grid(row=2, column=2, padx=5)
        ttk.Button(frame, text="Local PDB...", command=self.load_local_pdb).grid(row=2, column=3, padx=5)
        self.lbl_pdb_status = ttk.Label(frame, text=""); self.lbl_pdb_status.grid(row=3, column=3, sticky="w", padx=5)
        self.lbl_info = ttk.Label(frame, text="Protein not loaded", font=self.font_bold)
        self.lbl_info.grid(row=4, column=0, columnspan=14, sticky="w", pady=10)
        self.btn_clinvar = ttk.Button(frame, text="Analyze ClinVar", command=self.analyze_clinvar, state="disabled")
        self.btn_clinvar.grid(row=5, column=0, padx=5, pady=10, sticky="w")
        self.btn_save_txt = ttk.Button(frame, text="Save TXT", command=self.save_txt_report, state="disabled")
        self.btn_save_txt.grid(row=5, column=1, padx=5, pady=10, sticky="w")
        self.btn_save_html = ttk.Button(frame, text="Save HTML", command=self.save_html_report, state="disabled")
        self.btn_save_html.grid(row=5, column=2, padx=5, pady=10, sticky="w")
        self.btn_calibrate = ttk.Button(frame, text="Calibrate", command=self.calibrate_from_file)
        self.btn_calibrate.grid(row=5, column=3, padx=5, pady=10, sticky="w")
        ttk.Label(frame, text="Manual mutation:").grid(row=6, column=0, sticky="w", pady=5)
        self.entry_manual = tk.Entry(frame, width=12, bg='#3e3e3e', fg='#eff0f1', insertbackground='white', relief='solid', bd=1)
        self.entry_manual.grid(row=6, column=1, sticky="w", padx=10)
        self._add_context_menu(self.entry_manual)
        ttk.Button(frame, text="Check", command=self.analyze_manual).grid(row=6, column=2, padx=5)
        ttk.Label(frame, text="Second mutation:").grid(row=6, column=3, sticky="w", pady=5, padx=5)
        self.entry_manual2 = tk.Entry(frame, width=12, bg='#3e3e3e', fg='#eff0f1', insertbackground='white', relief='solid', bd=1)
        self.entry_manual2.grid(row=6, column=4, sticky="w", padx=5)
        self._add_context_menu(self.entry_manual2)
        ttk.Button(frame, text="Double", command=self.analyze_double).grid(row=6, column=5, padx=5)
        ttk.Button(frame, text="Open 3D in Browser", command=self.open_3d_browser).grid(row=6, column=6, padx=5)
        ttk.Label(frame, text="Batch mutations (one per line):").grid(row=7, column=0, sticky="w", pady=5)
        self.entry_batch = tk.Text(frame, width=40, height=6, bg='#3e3e3e', fg='#eff0f1', insertbackground='white', relief='solid', bd=1)
        self.entry_batch.grid(row=8, column=0, columnspan=2, sticky="w", padx=10, pady=5)
        self._add_context_menu(self.entry_batch)
        ttk.Button(frame, text="Batch Check", command=self.analyze_batch).grid(row=8, column=2, padx=5)
        ttk.Button(frame, text="Stats", command=self.show_analysis_stats).grid(row=8, column=3, padx=5)

        columns = ("pos","change","depth_old","depth_new","structure","classification",
                   "threshold","dG_bind","ddG","delta_hyd","delta_vol","ss_bonus","helix","strategy")
        self.tree = ttk.Treeview(frame, columns=columns, show="headings", height=18)
        for c, t in zip(columns, ["Pos.","Change","Old depth","New depth","Struct.","Class.",
                                  "Thresh. J","\u0394G bind.","\u0394\u0394G","\u0394Hyd","\u0394Vol","S-S","Helix","Strategy"]):
            self.tree.heading(c, text=t)
        for c, w in zip(columns, [40,70,80,80,70,80,70,70,70,50,50,40,40,140]):
            self.tree.column(c, width=w, anchor="center")
        self.tree.grid(row=9, column=0, columnspan=14, sticky="nsew", pady=10)
        frame.grid_rowconfigure(9, weight=1)
        self.tree_menu = tk.Menu(self.tree, tearoff=0, bg="#3e3e3e", fg="#eff0f1")
        self.tree_menu.add_command(label="Copy", command=self.copy_selected_rows)
        self.tree.bind("<Button-3>", self.show_tree_menu)
        self.tree.bind("<Control-c>", lambda e: self.copy_selected_rows())

    def build_graph_tab(self):
        self.graph_frame = ttk.Frame(self.tab_graph, padding=10); self.graph_frame.pack(fill="both", expand=True)

    def build_3d_tab(self):
        self.frame_3d = ttk.Frame(self.tab_3d, padding=10); self.frame_3d.pack(fill="both", expand=True)
        if not HAS_BIOPYTHON:
            ttk.Label(self.frame_3d, text="3D Structure unavailable.\nInstall biopython: pip install biopython",
                      foreground="#ff5555", background="#2d2d2d", font=self.font_default).pack(expand=True)
        else:
            ttk.Label(self.frame_3d, text="Load protein and PDB to see the 3D structure.", font=self.font_bold).pack()

    def build_stpc_tab(self):
        self.stpc_frame = ttk.Frame(self.tab_stpc, padding=10); self.stpc_frame.pack(fill="both", expand=True)
        ttk.Label(self.stpc_frame, text="STPC Crystal View \u2014 graph after analysis", font=self.font_bold).pack()

    def build_db_tab(self):
        frame = ttk.Frame(self.tab_db, padding=20); frame.pack(fill="both", expand=True)
        ttk.Label(frame, text="Local ClinVar database (5.7M mutations)", font=self.font_bold).pack(anchor="w", pady=5)
        ttk.Button(frame, text="Update ClinVar", command=self.update_db).pack(pady=15)
        self.lbl_db_status = ttk.Label(frame, text=""); self.lbl_db_status.pack()

    def build_reports_tab(self):
        frame = ttk.Frame(self.tab_reports, padding=20); frame.pack(fill="both", expand=True)
        ttk.Label(frame, text="Generated reports", font=self.font_bold).pack(anchor="w", pady=5)
        self.listbox_reports = tk.Listbox(frame, bg="#3e3e3e", fg="#eff0f1", height=20, font=self.font_default)
        self.listbox_reports.pack(fill="both", expand=True, pady=10)
        ttk.Button(frame, text="Open HTML", command=self.open_html_report).pack(pady=5)
        self.refresh_report_list()

    def build_about_tab(self):
        frame = ttk.Frame(self.tab_about, padding=30); frame.pack(fill="both", expand=True)
        tk.Label(frame, text=f"AllyMind {APP_VERSION}", font=self.font_title, fg="#3daee9", bg="#2d2d2d").pack(pady=10)
        tk.Label(frame, text=f"Author: {AUTHOR}\nIndependent researcher, Belarus\n\n{COPYRIGHT}",
                 justify="left", bg="#2d2d2d", fg="#eff0f1", font=self.font_default).pack()

    def build_isomer_tab(self):
        frame = ttk.Frame(self.tab_isomer, padding=20)
        frame.pack(fill="both", expand=True)
        ttk.Label(frame, text="Анализ стабильности ядерных изомеров", font=self.font_title).pack(pady=10)
        row1 = ttk.Frame(frame)
        row1.pack(fill="x", pady=5)
        ttk.Label(row1, text="A (массовое число):").pack(side="left")
        self.entry_A = ttk.Entry(row1, width=10)
        self.entry_A.pack(side="left", padx=5)
        ttk.Label(row1, text="Z (заряд):").pack(side="left", padx=10)
        self.entry_Z = ttk.Entry(row1, width=10)
        self.entry_Z.pack(side="left", padx=5)
        row2 = ttk.Frame(frame)
        row2.pack(fill="x", pady=5)
        ttk.Label(row2, text="Спин изомера (J):").pack(side="left")
        self.entry_spin_i = ttk.Entry(row2, width=10)
        self.entry_spin_i.pack(side="left", padx=5)
        ttk.Label(row2, text="Чётность изомера (+/-):").pack(side="left", padx=10)
        self.entry_parity_i = ttk.Entry(row2, width=5)
        self.entry_parity_i.pack(side="left", padx=5)
        row3 = ttk.Frame(frame)
        row3.pack(fill="x", pady=5)
        ttk.Label(row3, text="Спин осн. сост. (J):").pack(side="left")
        self.entry_spin_f = ttk.Entry(row3, width=10)
        self.entry_spin_f.pack(side="left", padx=5)
        ttk.Label(row3, text="Чётность осн. сост. (+/-):").pack(side="left", padx=10)
        self.entry_parity_f = ttk.Entry(row3, width=5)
        self.entry_parity_f.pack(side="left", padx=5)
        row4 = ttk.Frame(frame)
        row4.pack(fill="x", pady=5)
        ttk.Label(row4, text="Энергия перехода (МэВ):").pack(side="left")
        self.entry_deltaE = ttk.Entry(row4, width=10)
        self.entry_deltaE.pack(side="left", padx=5)
        btn_analyze = ttk.Button(frame, text="Анализировать изомер", command=self.analyze_isomer)
        btn_analyze.pack(pady=15)
        self.lbl_isomer_result = tk.Text(frame, height=4, bg="#3e3e3e", fg="#eff0f1", font=self.font_default, state="disabled", relief="solid", bd=1); self._add_context_menu(self.lbl_isomer_result)
        self.lbl_isomer_result.pack(pady=10, fill="x")
        ttk.Label(frame, text="Batch isomers (one per line: A Z J_i par_i J_f par_f dE_MeV):").pack(anchor="w", pady=(20,5))
        self.entry_batch_isomer = tk.Text(frame, width=60, height=6, bg='#3e3e3e', fg='#eff0f1', insertbackground='white', relief='solid', bd=1)
        self.entry_batch_isomer.pack(fill="x", pady=5)
        btn_batch_isomer = ttk.Button(frame, text="Batch Analyze Isomers", command=self.analyze_batch_isomer)
        btn_batch_isomer.pack(pady=5)
        self.lbl_batch_isomer_result = tk.Text(frame, height=10, bg="#3e3e3e", fg="#eff0f1", font=self.font_default, state="disabled", relief="solid", bd=1); self._add_context_menu(self.lbl_batch_isomer_result)
        self.lbl_batch_isomer_result.pack(pady=5, fill="x")

    def analyze_batch_isomer(self):
        text = self.entry_batch_isomer.get("1.0", "end-1c").strip()
        if not text:
            messagebox.showwarning("Batch Isomer", "Enter at least one isomer line."); return
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        results = []
        for ln in lines:
            parts = ln.split()
            if len(parts) < 6: continue
            try:
                A = int(parts[0]); Z = int(parts[1])
                J_i = float(parts[2]); par_i = parts[3]
                J_f = float(parts[4]); par_f = parts[5]
                dE_MeV = float(parts[6]) if len(parts) > 6 else 0.0
            except ValueError:
                results.append(f"Invalid line: {ln}"); continue
            dE_GeV = dE_MeV * 0.001
            L = abs(J_i - J_f)
            if L == 0: L = 1
            half_life = (1.0 / (dE_GeV ** (2*L + 1))) * 1e-12
            stability = "Stable" if half_life > 1e-9 else "Short-lived"
            results.append(f"A={A} Z={Z} J_i={J_i}{par_i} J_f={J_f}{par_f} dE={dE_MeV} MeV T1/2={half_life:.2e} s {stability}")
        self.lbl_batch_isomer_result.config(state="normal"); self.lbl_batch_isomer_result.delete("1.0", "end"); self.lbl_batch_isomer_result.insert("1.0", "\n".join(results)); self.lbl_batch_isomer_result.config(state="disabled")

    def analyze_isomer(self):
        try:
            A = int(self.entry_A.get())
            Z = int(self.entry_Z.get())
            spin_i = float(self.entry_spin_i.get())
            parity_i = self.entry_parity_i.get().strip()
            spin_f = float(self.entry_spin_f.get())
            parity_f = self.entry_parity_f.get().strip()
            delta_E_MeV = float(self.entry_deltaE.get())
        except ValueError:
            messagebox.showerror("Ошибка", "Проверьте правильность введённых чисел")
            return
        delta_E_GeV = delta_E_MeV * 0.001
        L = abs(spin_i - spin_f)
        if L == 0: L = 1
        half_life = (1.0 / (delta_E_GeV ** (2*L + 1))) * 1e-12
        text = f"Период полураспада: {half_life:.2e} с\n"
        text += f"Энергия перехода: {delta_E_MeV:.2f} МэВ\n"
        text += f"Мультипольность L: {L}\n"
        text += f"Класс: {'Stable' if half_life > 1e-9 else 'Short-lived'}"
        self.lbl_isomer_result.config(state="normal"); self.lbl_isomer_result.delete("1.0", "end"); self.lbl_isomer_result.insert("1.0", text); self.lbl_isomer_result.config(state="disabled")

    def load_protein(self):
        query = self.entry_id.get().strip().upper()
        if not query:
            messagebox.showwarning("Input", "Enter UniProt ID or gene symbol"); return
        self.set_status("Loading...")
        try:
            seq, uid = allymind_core.load_sequence(query)
            if seq:
                self.uniprot_id = uid
                self.gene_symbol = query
                self.sequence = seq
                self.ss_map = None
                self.pdb_path = None
                self.metal_sites = []
                self.metal_coordinated_residues = set()
                self.dimer_interface_residues = set()
                self.lbl_pdb_status.config(text="")
                self.entry_pdb.delete(0, tk.END)
                self._update_3d_structure()
                self._on_loaded()
            else:
                messagebox.showerror("Error", f"Could not find protein '{query}'.")
        except Exception as e:
            messagebox.showerror("Loading Error", str(e))
        self.set_status("Ready")

    def load_fasta(self):
        path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa")])
        if not path: return
        try:
            with open(path, encoding='utf-8') as f:
                self.sequence = ''.join(line.strip() for line in f if not line.startswith('>'))
            self.entry_pdb.delete(0, tk.END)
        except Exception as e:
            messagebox.showerror("Error", str(e)); return
        gene = self.entry_id.get().strip().upper()
        if not gene:
            messagebox.showwarning("Input", "Enter gene symbol"); return
        self.gene_symbol = gene
        self.uniprot_id = gene
        self.ss_map = None
        self.pdb_path = None
        self.metal_sites = []
        self.metal_coordinated_residues = set()
        self.dimer_interface_residues = set()
        self.lbl_pdb_status.config(text="")
        self._update_3d_structure()
        self._on_loaded()

    def load_pdb(self):
        self.pdb_path = None
        self.ss_map = None
        raw_pdb = self.entry_pdb.get().strip().upper()
        if not raw_pdb:
            messagebox.showwarning("Input", "Enter PDB ID (e.g., 3D09)")
            return
        if not self.sequence:
            messagebox.showwarning("PDB", "Load a protein first.")
            return
        if '_' in raw_pdb:
            parts = raw_pdb.split('_')
            pdb_id = parts[0]
            target_chain = parts[1]
        else:
            pdb_id = raw_pdb
            target_chain = 'A'
        self.set_status(f"Loading PDB {pdb_id} (chain {target_chain})...")
        try:
            self.pdb_path = pdb_parser.download_pdb(pdb_id)
            if not self.pdb_path:
                raise Exception("Could not download PDB file.")
            self.ss_map = pdb_parser.extract_aligned_ss_map(self.sequence, self.pdb_path, chain=target_chain)
            self.metal_sites = pdb_parser.extract_metal_sites(self.pdb_path, target_chain=target_chain)
            self.metal_coordinated_residues = cascade.detect_metal_coordinations(self.metal_sites, self.pdb_path)
            self.dimer_interface_residues = cascade.detect_dimer_interface(self.pdb_path)
            self.lbl_pdb_status.config(text=f"OK {pdb_id} chain {target_chain} loaded, metals: {len(self.metal_sites)}")
            self.set_status(f"PDB {pdb_id} loaded. Metal sites: {len(self.metal_sites)}")
            self._update_3d_structure()
        except Exception as e:
            messagebox.showerror("PDB Error", str(e))
            self.ss_map = None
            self.pdb_path = None
            self.metal_sites = []
            self.metal_coordinated_residues = set()
            self.dimer_interface_residues = set()
            self.lbl_pdb_status.config(text="")
            self.set_status("PDB loading error")

    def load_local_pdb(self):
        self.pdb_path = None; self.ss_map = None
        path = filedialog.askopenfilename(filetypes=[("PDB files","*.pdb"),("All files","*.*")])
        if not path: return
        if not self.sequence:
            messagebox.showwarning("PDB", "Load a protein first."); return
        self.set_status(f"Loading local PDB {os.path.basename(path)}...")
        try:
            self.pdb_path = path
            self.ss_map = pdb_parser.extract_aligned_ss_map(self.sequence, path)
            self.metal_sites = pdb_parser.extract_metal_sites(path)
            self.metal_coordinated_residues = cascade.detect_metal_coordinations(self.metal_sites, path)
            self.dimer_interface_residues = cascade.detect_dimer_interface(path)
            self.lbl_pdb_status.config(text=f"OK {os.path.basename(path)} loaded, metals: {len(self.metal_sites)}")
            self.set_status(f"Local PDB loaded. Metal sites: {len(self.metal_sites)}")
            self._update_3d_structure()
        except Exception as e:
            messagebox.showerror("PDB Error", str(e))
            self.ss_map = None
            self.pdb_path = None
            self.metal_sites = []
            self.metal_coordinated_residues = set()
            self.dimer_interface_residues = set()
            self.lbl_pdb_status.config(text="")
            self.set_status("PDB loading error")

    def analyze_manual(self):
        start_time = time.time()
        mut = self.entry_manual.get().strip()
        if not mut or len(mut) < 3: return
        wt, pos, alt = mut[0].upper(), int(mut[1:-1])-1, mut[-1].upper()
        if not self.sequence: return
        report = allymind_core.analyze_mutation(self.sequence, mut.upper(), ss_map=self.ss_map)
        trap_old, depth_old = check_trap_old(self.sequence, pos, self.ss_map)
        depth_new = report['depth']
        if wt == alt:
            depth_new = depth_old
        exposed = cascade.is_trap_exposed(self.sequence, pos, wt, alt, self.ss_map)
        classification = "open" if exposed else "hidden"
        thr = report.get('therapeutic_threshold', 0) or 0
        dG = (report.get('kinetics') or {}).get('delta_G', 0) or 0
        ddG = report.get('delta_delta_G', 0) or 0
        ss = report.get('secondary_structure', '\u2014') or '\u2014'
        comp = report.get('chemical_components', {})

        strategy = determine_strategy(
            mut.upper(), depth_old, depth_new, classification, thr,
            self.uniprot_id, self.entry_pdb.get().strip().upper(),
            self.sequence, self.ss_map, self.pdb_path,
            self.metal_coordinated_residues,
            dimer_interface_residues=self.dimer_interface_residues
        )

        self.tree.delete(*self.tree.get_children())
        self.tree.insert("", "end", values=(
            pos+1, f"{wt}\u2192{alt}", f"{depth_old:.3f}", f"{depth_new:.3f}",
            ss, classification, f"{thr:.4f}", f"{dG:.4f}", f"{ddG:.4f}",
            f"{comp.get('delta_hyd', 0):.4f}",
            f"{comp.get('delta_vol', 0):.4f}",
            f"{comp.get('ss_bonus', 0):.4f}",
            f"{comp.get('helix_penalty', 0):.4f}",
            strategy
        ))
        self.current_mutation = mut
        self.plot_landscape_for_mutation(pos)
        self._update_3d_structure()
        self._plot_stpc_graph()
        self._increment_counter()
        elapsed = (time.time() - start_time) * 1000
        self.set_status(f"Old: {'trap' if trap_old else 'none'}, New: {classification} | {elapsed:.1f} ms")

    def analyze_clinvar(self):
        start_time = time.time()
        if not self.sequence or not self.gene_symbol: return
        self.set_status("Analyzing ClinVar...")
        try:
            conn = sqlite3.connect(allymind_core.DB_CLINVAR)
            rows = conn.execute("""SELECT protein_position, ref_aa, alt_aa FROM mutations
                                   WHERE gene=? AND protein_position IS NOT NULL
                                   AND clinical_significance LIKE '%athogenic%'""", (self.gene_symbol,)).fetchall()
            conn.close()
        except Exception as e:
            messagebox.showerror("DB Error", str(e)); return
        mutations = [(p, r, a) for p, r, a in rows if p is not None]
        if not mutations:
            messagebox.showinfo("ClinVar", f"No mutations found for {self.gene_symbol}"); return
        self.show_results(mutations)
        self._increment_counter()
        elapsed = (time.time() - start_time) * 1000
        self.set_status(f"ClinVar done: {len(mutations)} mutations | {elapsed:.1f} ms")

    def show_results(self, mutations):
        self.tree.delete(*self.tree.get_children())
        graph_data = []
        traps = 0
        for pos, ref, alt in mutations:
            trap_old, depth_old = check_trap_old(self.sequence, pos, self.ss_map)
            report = allymind_core.analyze_mutation(self.sequence, f"{ref}{pos}{alt}", ss_map=self.ss_map)
            trap_new, depth_new = report['trap'], report['depth']
            if ref == alt: depth_new = depth_old
            exposed = cascade.is_trap_exposed(self.sequence, pos, ref, alt, self.ss_map)
            if trap_new: traps += 1
            classification = "open" if exposed else "hidden"
            thr = report.get('therapeutic_threshold', 0) or 0
            dG = (report.get('kinetics') or {}).get('delta_G', 0) or 0
            ddG = report.get('delta_delta_G', 0) or 0
            ss = report.get('secondary_structure', '\u2014') or '\u2014'
            comp = report.get('chemical_components', {})
            strategy = determine_strategy(
                f"{ref}{pos}{alt}", depth_old, depth_new, classification, thr,
                self.uniprot_id, self.entry_pdb.get().strip().upper(),
                self.sequence, self.ss_map, self.pdb_path,
                self.metal_coordinated_residues,
                dimer_interface_residues=self.dimer_interface_residues
            )
            self.tree.insert("", "end", values=(
                pos+1, f"{ref}\u2192{alt}", f"{depth_old:.3f}", f"{depth_new:.3f}", ss, classification,
                f"{thr:.4f}", f"{dG:.4f}", f"{ddG:.4f}",
                f"{comp.get('delta_hyd', 0):.4f}", f"{comp.get('delta_vol', 0):.4f}",
                f"{comp.get('ss_bonus', 0):.4f}", f"{comp.get('helix_penalty', 0):.4f}", strategy
            ))
            graph_data.append((pos+1, depth_old or 0, depth_new, exposed, ref, alt))
        for w in self.graph_frame.winfo_children(): w.destroy()
        if graph_data: self.plot_distribution_in_frame(graph_data)
        self._plot_stpc_graph()

    def analyze_batch(self):
        start_time = time.time()
        if not self.sequence:
            messagebox.showwarning("Batch", "Load a protein first."); return
        text = self.entry_batch.get("1.0", "end-1c").strip()
        if not text:
            messagebox.showwarning("Batch", "Enter at least one mutation."); return
        lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
        for line in lines:
            parts = line.split()
            if len(parts) >= 3:
                gene, pdb_id, mut = parts[0].upper(), parts[1].lower(), parts[2].upper()
            elif len(parts) == 1:
                mut = parts[0].upper()
                gene = self.gene_symbol
                pdb_id = self.entry_pdb.get().strip().lower()
            else: continue
            if gene != self.gene_symbol or (pdb_id and pdb_id != self.entry_pdb.get().strip().lower()):
                try:
                    seq, uid = allymind_core.load_sequence(gene)
                    if not seq: continue
                    self.gene_symbol = gene; self.uniprot_id = uid; self.sequence = seq
                    self.entry_id.delete(0, tk.END); self.entry_id.insert(0, gene)
                    self.ss_map = None; self.pdb_path = None
                    self.metal_sites = []; self.metal_coordinated_residues = set(); self.dimer_interface_residues = set()
                    if pdb_id:
                        pdb_file = pdb_parser.download_pdb(pdb_id)
                        if pdb_file:
                            self.pdb_path = pdb_file
                            self.ss_map = pdb_parser.extract_aligned_ss_map(self.sequence, pdb_file)
                            self.metal_sites = pdb_parser.extract_metal_sites(pdb_file)
                            self.metal_coordinated_residues = cascade.detect_metal_coordinations(self.metal_sites, pdb_file)
                            self.dimer_interface_residues = cascade.detect_dimer_interface(pdb_file)
                            self.entry_pdb.delete(0, tk.END); self.entry_pdb.insert(0, pdb_id.upper())
                except: continue
            if len(mut) < 3: continue
            wt, alt = mut[0], mut[-1]
            try: pos = int(mut[1:-1]) - 1
            except: continue
            try:
                report = allymind_core.analyze_mutation(self.sequence, mut, ss_map=self.ss_map)
                trap_old, depth_old = check_trap_old(self.sequence, pos, self.ss_map)
                depth_new = report['depth']
                if wt == alt: depth_new = depth_old
                exposed = cascade.is_trap_exposed(self.sequence, pos, wt, alt, self.ss_map)
                classification = "open" if exposed else "hidden"
                thr = report.get('therapeutic_threshold', 0) or 0
                dG = (report.get('kinetics') or {}).get('delta_G', 0) or 0
                ddG = report.get('delta_delta_G', 0) or 0
                ss = report.get('secondary_structure', '\u2014') or '\u2014'
                comp = report.get('chemical_components', {})
                strategy = determine_strategy(
                    mut, depth_old, depth_new, classification, thr,
                    self.uniprot_id, pdb_id.upper() if pdb_id else self.entry_pdb.get().strip().upper(),
                    self.sequence, self.ss_map, self.pdb_path,
                    self.metal_coordinated_residues,
                    dimer_interface_residues=self.dimer_interface_residues
                )
                self.tree.insert("", "end", values=(
                    pos+1, f"{wt}\u2192{alt}", f"{depth_old:.3f}", f"{depth_new:.3f}", ss, classification,
                    f"{thr:.4f}", f"{dG:.4f}", f"{ddG:.4f}",
                    f"{comp.get('delta_hyd', 0):.4f}", f"{comp.get('delta_vol', 0):.4f}",
                    f"{comp.get('ss_bonus', 0):.4f}", f"{comp.get('helix_penalty', 0):.4f}", strategy
                ))
            except Exception as e:
                print(f"Batch error for {mut}: {e}")
        self._increment_counter()
        elapsed = (time.time() - start_time) * 1000
        self.set_status(f"Batch done: {len(lines)} lines | {elapsed:.1f} ms")
        self._plot_stpc_graph()
        self._update_3d_structure()

    def _on_loaded(self):
        self.lbl_info.config(text=f"{self.gene_symbol} (UniProt: {self.uniprot_id}), length {len(self.sequence)} a.a.")
        self.btn_clinvar.config(state="normal")
        self.btn_save_txt.config(state="normal")
        self.btn_save_html.config(state="normal")
        self._update_3d_structure()
        self.set_status("Ready")

    def plot_distribution_in_frame(self, data):
        pos = [d[0] for d in data]; old = [d[1] for d in data]; new = [d[2] for d in data]
        exp = [d[3] for d in data]; codes = [f"{d[4]}{d[0]}{d[5]}" for d in data]
        fig, ax = plt.subplots(figsize=(7,4))
        ax.scatter(pos, old, c='red', marker='o', alpha=0.5, label='Old model')
        ax.scatter(pos, new, c='#3daee9', marker='s', alpha=0.5, label='New model')
        open_pos = [p for p,e in zip(pos,exp) if e]; open_depths = [d for d,e in zip(new,exp) if e]
        if open_pos: ax.scatter(open_pos, open_depths, c='none', edgecolors='gold', s=80, linewidths=1.5, label='Open traps')
        ax.set_xlabel('Position'); ax.set_ylabel('Depth (HB)'); ax.legend(); ax.grid(alpha=0.3)
        fig.tight_layout()
        annot = ax.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
                            bbox=dict(boxstyle="round",fc="w"), arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        def hover(event):
            if event.inaxes != ax: annot.set_visible(False); fig.canvas.draw_idle(); return
            min_dist, closest = float('inf'), None
            for i,(x,y) in enumerate(zip(pos,new)):
                d = (event.xdata-x)**2 + (event.ydata-y)**2
                if d < min_dist: min_dist = d; closest = i
            if closest is not None and min_dist<100:
                annot.xy = (pos[closest], new[closest])
                annot.set_text(f"{codes[closest]}\nOld: {old[closest]:.3f} HB\nNew: {new[closest]:.3f} HB")
                annot.set_visible(True)
            else: annot.set_visible(False)
            fig.canvas.draw_idle()
        fig.canvas.mpl_connect("motion_notify_event", hover)
        canvas = FigureCanvasTkAgg(fig, master=self.graph_frame); canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self.graph_frame); toolbar.update(); toolbar.pack(side='bottom', fill='x')
        self._add_zoom_and_pan(fig, ax)
        self._last_fig = fig; self._last_scatter_data = data

    def plot_landscape_for_mutation(self, pos):
        for w in self.graph_frame.winfo_children(): w.destroy()
        if not self.sequence or not self.current_mutation: return
        J_old, _ = cascade.secondary_structure_J(self.sequence, pos, self.ss_map)
        J_new, dom_rad = cascade.secondary_structure_J(self.sequence, pos, self.ss_map)
        start = max(0, pos-dom_rad); end = min(len(self.sequence), pos+dom_rad+1)
        dom_len = end-start; healthy = len(self.sequence)-dom_len
        n_vals = np.arange(0, len(self.sequence)+1)
        F_wt = np.where(n_vals <= healthy,
                        cascade.monomer_free_energy(n_vals, 0, cascade.J_NAT, len(self.sequence)),
                        cascade.monomer_free_energy(healthy, n_vals-healthy, cascade.J_NAT, len(self.sequence)))
        F_old = np.where(n_vals <= healthy,
                         cascade.monomer_free_energy(n_vals, 0, J_old, len(self.sequence)),
                         cascade.monomer_free_energy(healthy, n_vals-healthy, J_old, len(self.sequence)))
        F_new = np.where(n_vals <= healthy,
                         cascade.monomer_free_energy(n_vals, 0, J_new, len(self.sequence)),
                         cascade.monomer_free_energy(healthy, n_vals-healthy, J_new, len(self.sequence)))
        fig, ax = plt.subplots(figsize=(6,3.5))
        ax.plot(n_vals/len(self.sequence), F_wt, '#eff0f1', lw=2, label='WT')
        ax.plot(n_vals/len(self.sequence), F_old, 'r--', lw=2, label=f'Old (J={J_old:.2f})')
        ax.plot(n_vals/len(self.sequence), F_new, '#3daee9', lw=2, label=f'New (J={J_new:.2f})')
        ax.set_xlabel('Folded fraction'); ax.set_ylabel('Free energy (HB)'); ax.legend(); ax.grid(alpha=0.3)
        fig.tight_layout()
        canvas = FigureCanvasTkAgg(fig, master=self.graph_frame); canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self.graph_frame); toolbar.update(); toolbar.pack(side='bottom', fill='x')
        self._add_zoom_and_pan(fig, ax)

    def save_txt_report(self):
        if not self.tree.get_children(): return
        fn = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files","*.txt")])
        if not fn: return
        with open(fn, 'w', encoding='utf-8') as f:
            f.write(f"AllyMind {APP_VERSION} Report: {self.gene_symbol}\n{COPYRIGHT}\n")
            f.write(f"Length: {len(self.sequence)} a.a. PDB: {self.entry_pdb.get() or 'none'}\n")
            hdr = f"{'Pos.':<6} {'Change':<10} {'Old depth':<10} {'New depth':<10} {'Struct.':<10} {'Class':<10} {'Thresh. J':<10} {'ΔG':<10} {'ΔΔG':<10} {'ΔHyd':<10} {'ΔVol':<10} {'S-S':<10} {'Helix':<10} {'Strategy':<20}"
            f.write(hdr + "\n")
            for child in self.tree.get_children():
                vals = self.tree.item(child, 'values')
                f.write(f"{vals[0]:<6} {vals[1]:<10} {vals[2]:<10} {vals[3]:<10} {vals[4]:<10} {vals[5]:<10} {vals[6]:<10} {vals[7]:<10} {vals[8]:<10} {vals[9]:<10} {vals[10]:<10} {vals[11]:<10} {vals[12]:<10} {vals[13]:<20}\n")
        messagebox.showinfo("TXT Report", f"Saved to {fn}")

    def save_html_report(self):
        html = self._build_html_content()
        if not html: return
        os.makedirs("reports", exist_ok=True)
        fn = f"reports/report_{self.gene_symbol}_{datetime.datetime.now():%Y%m%d_%H%M%S}.html"
        with open(fn, 'w', encoding='utf-8') as f: f.write(html)
        messagebox.showinfo("HTML Report", f"Saved to {fn}"); self.refresh_report_list()

    def _build_html_content(self):
        if not self.tree.get_children(): return None
        rows = [self.tree.item(i,'values') for i in self.tree.get_children()]
        total = len(rows); open_traps = sum(1 for r in rows if 'open' in str(r[5]))
        chart_data = []
        if self._last_scatter_data:
            for p,o,n,e,ref,alt in self._last_scatter_data:
                chart_data.append({'x':p,'y_old':o,'y_new':n,'exposed':e,'code':f"{ref}{p}{alt}"})
        pdb_id = self.entry_pdb.get().strip().upper() or "1UBQ"
        html = f"""<!DOCTYPE html><html lang="en">
<head><meta charset="UTF-8"><title>Report: {self.gene_symbol}</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0"></script>
<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-zoom@2.0.1"></script>
<style>
body{{font-family:'Segoe UI',sans-serif;background:#1a1a1a;color:#eff0f1;padding:20px}}
h1{{text-align:center;border-bottom:2px solid #3daee9;padding-bottom:10px;color:#3daee9}}
.summary{{display:flex;justify-content:space-around;margin:20px 0}}
.summary div{{background:#2d2d2d;padding:15px;border-radius:8px;text-align:center}}
table{{width:100%;border-collapse:collapse;margin-top:20px}}
th,td{{border:1px solid #3daee9;padding:8px;font-size:0.9rem}}
th{{background:#333}}.open{{color:#ff5555;font-weight:bold}}.hidden{{color:#55ff55}}
.chart-container{{width:800px;height:400px;margin:20px auto}}
.chart-controls{{text-align:center;margin-bottom:20px}}
.chart-controls button{{background:#3daee9;color:white;border:none;padding:8px 16px;margin:0 5px;border-radius:4px;cursor:pointer}}
.chart-controls button:hover{{background:#2980b9}}
.footer{{text-align:center;margin-top:40px;font-size:0.8rem;color:#555}}
</style></head>
<body>
<h1>AllyMind {APP_VERSION} \u2014 Report for {self.gene_symbol}</h1>
<p style="text-align:center">Date: {datetime.datetime.now():%d.%m.%Y %H:%M} | Length: {len(self.sequence)} a.a. | PDB: {pdb_id}</p>
<div class="summary"><div><b>Total</b><br>{total}</div><div><b>Open</b><br>{open_traps}</div><div><b>Hidden</b><br>{total-open_traps}</div></div>
<div class="chart-container"><canvas id="trapChart"></canvas></div>
<div class="chart-controls">
<button onclick="chart.zoom(1.2)">🔍 +</button>
<button onclick="chart.zoom(0.8)">🔍 −</button>
<button onclick="chart.resetZoom()">↺ Reset</button>
</div>
<script>
const rawData = {json.dumps(chart_data)};
const ctx = document.getElementById('trapChart').getContext('2d');
window.chart = new Chart(ctx,{{type:'scatter',data:{{datasets:[
{{label:'Old model',data:rawData.map(d=>({{x:d.x,y:d.y_old}})),backgroundColor:'rgba(255,0,0,0.5)',pointRadius:4}},
{{label:'New model',data:rawData.map(d=>({{x:d.x,y:d.y_new}})),backgroundColor:'rgba(61,174,233,0.5)',pointRadius:4}},
{{label:'Open traps',data:rawData.filter(d=>d.exposed).map(d=>({{x:d.x,y:d.y_new}})),backgroundColor:'gold',borderColor:'orange',pointRadius:6}}]}},
options:{{responsive:true,maintainAspectRatio:false,
plugins:{{zoom:{{zoom:{{wheel:{{enabled:true}},pinch:{{enabled:true}},drag:{{enabled:true}}}},pan:{{enabled:true}}}}}},
scales:{{x:{{title:{{display:true,text:'Position'}}}},y:{{title:{{display:true,text:'Depth (HB)'}}}}}}}}}});
</script>
<h2>Mutation table</h2>
<table><tr><th>Pos.</th><th>Change</th><th>Old depth</th><th>New depth</th><th>Struct.</th><th>Class.</th><th>Thresh. J</th><th>ΔG bind.</th><th>ΔΔG</th><th>ΔHyd</th><th>ΔVol</th><th>S‑S</th><th>Helix</th><th>Strategy</th></tr>"""
        for r in rows:
            cls = 'open' if 'open' in str(r[5]) else 'hidden'
            html += f"<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td>{r[3]}</td><td>{r[4]}</td><td class='{cls}'>{r[5]}</td><td>{r[6]}</td><td>{r[7]}</td><td>{r[8]}</td><td>{r[9]}</td><td>{r[10]}</td><td>{r[11]}</td><td>{r[12]}</td><td>{r[13]}</td></tr>\n"
        html += "</table><div class='footer'><p>AllyMind 10.2 FINAL \u00a9 Valentin Lebedkin, 2026. Confidential.</p></div></body></html>"
        return html

    def open_3d_browser(self):
        if not self.current_mutation:
            messagebox.showwarning("3D", "Сначала выполните анализ мутации (Check)."); return
        pdb_id = self.entry_pdb.get().strip().upper() or "2C9V"
        try: pos = int(self.current_mutation[1:-1])
        except: pos = 4
        html = f"""<!DOCTYPE html><html><head><meta charset="UTF-8">
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>body,html{{margin:0;padding:0;width:100%;height:100%;overflow:hidden;background:#1e1e24;}}#viewer{{width:100%;height:100%;position:relative;}}</style>
</head><body><div id="viewer"></div><script>
let viewer = $3Dmol.createViewer('viewer', {{ backgroundColor:'#1e1e24' }});
$3Dmol.download('pdb:{pdb_id}', viewer, {{}}, function() {{
    viewer.setStyle({{}}, {{ cartoon: {{ color:'spectrum' }} }});
    viewer.setStyle({{ chain:'A', resi:{pos} }}, {{ sphere: {{ color:'#ff3344', radius:1.2 }}, cartoon: {{ color:'#ff3344' }} }});
    viewer.addLabel("A:{pos}", {{ position:{{ chain:'A', resi:{pos} }}, backgroundColor:'#222', fontColor:'white', fontSize:14 }});
    viewer.zoomTo(); viewer.render();
}});</script></body></html>"""
        path = os.path.join(tempfile.gettempdir(), f"protein_{pdb_id}_{pos}.html")
        with open(path, 'w', encoding='utf-8') as f: f.write(html)
        webbrowser.open('file:///' + path.replace('\\', '/'))

    def update_db(self):
        self.set_status("Updating ClinVar..."); subprocess.Popen([sys.executable, "clinvar_importer_v2.py"])
        self.lbl_db_status.config(text="Запущено обновление (см. консоль)")

    def refresh_report_list(self):
        self.listbox_reports.delete(0, tk.END)
        if os.path.exists("reports"):
            for f in os.listdir("reports"):
                if f.endswith(".html"): self.listbox_reports.insert(tk.END, os.path.join("reports", f))

    def open_html_report(self):
        sel = self.listbox_reports.curselection()
        if sel: os.startfile(self.listbox_reports.get(sel[0]))

    def calibrate_from_file(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files","*.csv"),("Text files","*.txt")])
        if not path: return
        try:
            import csv
            exp, calc = [], []
            with open(path,'r') as fh:
                reader = csv.reader(fh)
                next(reader,None)
                for row in reader:
                    if len(row)>=3:
                        mut = row[0].strip(); exp_ddg = float(row[1]); hb = float(row[2])
                        exp.append(exp_ddg); calc.append(hb)
            result = cascade.calibrate_to_protherm(exp, calc)
            if result:
                self.calibration_scale, self.calibration_intercept = result
                messagebox.showinfo("Калибровка", f"Масштаб: {self.calibration_scale:.4f} ккал/моль на HB\nСдвиг: {self.calibration_intercept:.4f}")
            else: messagebox.showerror("Ошибка", "Недостаточно данных для калибровки.")
        except Exception as e:
            messagebox.showerror("Ошибка", str(e))

    def set_status(self, text): self.status.config(text=text); self.root.update_idletasks()

    def _update_3d_structure(self):
        for w in self.frame_3d.winfo_children(): w.destroy()
        if not HAS_BIOPYTHON:
            ttk.Label(self.frame_3d, text="3D Structure unavailable.\nInstall biopython: pip install biopython",
                      foreground="#ff5555", background="#2d2d2d", font=self.font_default).pack(expand=True)
            return
        if not self.sequence or not self.ss_map or self.pdb_path is None:
            ttk.Label(self.frame_3d, text="Load a protein and PDB structure first.",
                      foreground="#888888", background="#2d2d2d").pack(expand=True)
            return
        pdb_file = self.pdb_path
        if not pdb_file:
            pdb_id = self.entry_pdb.get().strip().upper()
            pdb_file = pdb_parser.download_pdb(pdb_id)
        if not pdb_file:
            ttk.Label(self.frame_3d, text="PDB file not available.", foreground="#ff5555").pack(expand=True)
            return
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('prot', pdb_file)
        except Exception as e:
            ttk.Label(self.frame_3d, text=f"Error parsing PDB: {e}", foreground="#ff5555").pack(expand=True)
            return
        ca = []
        for model in structure:
            for chain in model:
                for res in chain:
                    if 'CA' in res: ca.append(res['CA'].coord)
        if not ca:
            ttk.Label(self.frame_3d, text="No CA atoms found in PDB.", foreground="#ff5555").pack(expand=True)
            return
        ca = np.array(ca)
        cf = cascade.predict_ss_chou_fasman(self.sequence)
        colors = []
        for i in range(len(ca)):
            if self.ss_map and i < len(self.ss_map) and self.ss_map[i] in ('H','E','C'):
                ss = self.ss_map[i]
            elif i < len(cf):
                ss = cf[i]
            else:
                ss = 'C'
            colors.append('#ff5555' if ss == 'H' else '#55ff55' if ss == 'E' else '#888888')
        fig = plt.figure(figsize=(7,6), facecolor='#2d2d2d')
        ax = fig.add_subplot(111, projection='3d')
        ax.set_facecolor('#2d2d2d')
        ax.set_title(f'3D Structure (PDB: {self.entry_pdb.get() or "local"})', color='white')
        ax.scatter(ca[:,0], ca[:,1], ca[:,2], c=colors, s=10, alpha=0.9, edgecolors='none')
        for i in range(len(ca)-1):
            ax.plot(ca[i:i+2,0], ca[i:i+2,1], ca[i:i+2,2], color='#555555', linewidth=0.5, antialiased=True)
        if self.current_mutation:
            try:
                pos = int(self.current_mutation[1:-1])-1
                if 0 <= pos < len(ca):
                    ax.scatter([ca[pos,0]], [ca[pos,1]], [ca[pos,2]], c='gold', s=100, edgecolors='black', linewidth=2)
                    ax.text(ca[pos,0], ca[pos,1], ca[pos,2]+1.5, self.current_mutation, color='white', fontsize=8, fontweight='bold', bbox=dict(facecolor='black', alpha=0.5, edgecolor='none'))
            except: pass
        ax.grid(False); ax.axis('off')
        canvas = FigureCanvasTkAgg(fig, master=self.frame_3d); canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self.frame_3d); toolbar.update(); toolbar.pack(side='bottom', fill='x')
        self._add_3d_navigation_rot(fig, ax)

    def _plot_stpc_graph(self):
        for w in self.stpc_frame.winfo_children():
            if w != self.stpc_frame.winfo_children()[0]: w.destroy()
        if not self.sequence: return
        N = len(self.sequence)
        coords = None
        if self.pdb_path and HAS_BIOPYTHON:
            try:
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure('prot', self.pdb_path)
                ca = []
                for model in structure:
                    for chain in model:
                        for res in chain:
                            if 'CA' in res: ca.append(res['CA'].coord)
                if ca and len(ca) >= N: coords = np.array(ca[:N])
            except: pass
        if coords is None:
            t = np.linspace(0, 4*np.pi, N)
            coords = np.column_stack((np.sin(t)*20, np.cos(t)*20, np.linspace(0, N*0.3, N)))
        fig = plt.figure(figsize=(7,6), facecolor='#2d2d2d')
        ax = fig.add_subplot(111, projection='3d')
        ax.set_facecolor('#2d2d2d')
        for i in range(N-1):
            ax.plot(coords[i:i+2,0], coords[i:i+2,1], coords[i:i+2,2], color='#555555', linewidth=0.5, antialiased=True)
        danger = []
        for i in range(N):
            if self.ss_map and i < len(self.ss_map):
                J, _ = cascade.secondary_structure_J(self.sequence, i, self.ss_map)
                report = allymind_core.analyze_mutation(self.sequence, f"A{i+1}V", ss_map=self.ss_map)
                depth = report.get('depth',0); thr = report.get('therapeutic_threshold',0)
                danger.append('#ff4444' if depth > 3.0 else '#ffaa00' if depth > thr else '#888888')
            else: danger.append('#888888')
        ax.scatter(coords[:,0], coords[:,1], coords[:,2], c=danger, s=15, alpha=0.9)
        for i in range(0, N, 10):
            ax.text(coords[i,0], coords[i,1], coords[i,2], str(i+1), color='white', fontsize=6)
        if self.current_mutation:
            try:
                pos = int(self.current_mutation[1:-1])-1
                if 0 <= pos < N:
                    report = allymind_core.analyze_mutation(self.sequence, self.current_mutation, ss_map=self.ss_map)
                    J = report.get('J_local',0); depth = report.get('depth',0)
                    ax.text(coords[pos,0], coords[pos,1], coords[pos,2]+1.5, f'J={J:.3f}\nd={depth:.3f}', color='white', fontsize=8, fontweight='bold', bbox=dict(facecolor='black', alpha=0.5, edgecolor='none'))
            except: pass
        ax.set_title('STPC Crystal Graph', color='white'); ax.grid(False); ax.axis('off')
        canvas = FigureCanvasTkAgg(fig, master=self.stpc_frame); canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(canvas, self.stpc_frame); toolbar.update(); toolbar.pack(side='bottom', fill='x')
        self._add_3d_navigation_rot(fig, ax)

    def _add_zoom_and_pan(self, fig, ax):
        def zoom_fun(event):
            if event.inaxes != ax: return
            cur_xlim, cur_ylim = ax.get_xlim(), ax.get_ylim()
            xdata, ydata = event.xdata, event.ydata
            if xdata is None or ydata is None: return
            scale = 1.2 if event.step>0 else 0.8
            nw = (cur_xlim[1]-cur_xlim[0])*scale; nh = (cur_ylim[1]-cur_ylim[0])*scale
            rx = (cur_xlim[1]-xdata)/(cur_xlim[1]-cur_xlim[0]); ry = (cur_ylim[1]-ydata)/(cur_ylim[1]-cur_ylim[0])
            ax.set_xlim([xdata-nw*(1-rx), xdata+nw*rx]); ax.set_ylim([ydata-nh*(1-ry), ydata+nh*ry])
            fig.canvas.draw_idle()
        def pan_press(event):
            if event.inaxes != ax: return
            self._pan_active = True; self._pan_start_x = event.xdata; self._pan_start_y = event.ydata
        def pan_release(event): self._pan_active = False
        def pan_motion(event):
            if not self._pan_active or event.inaxes != ax: return
            if event.xdata is None or event.ydata is None: return
            dx = self._pan_start_x - event.xdata; dy = self._pan_start_y - event.ydata
            ax.set_xlim(ax.get_xlim()+dx); ax.set_ylim(ax.get_ylim()+dy)
            fig.canvas.draw_idle()
        fig.canvas.mpl_connect('scroll_event', zoom_fun)
        fig.canvas.mpl_connect('button_press_event', pan_press)
        fig.canvas.mpl_connect('button_release_event', pan_release)
        fig.canvas.mpl_connect('motion_notify_event', pan_motion)

    def _add_3d_navigation_rot(self, fig, ax):
        def zoom_fun(event):
            if event.inaxes != ax: return
            cur_xlim, cur_ylim = ax.get_xlim(), ax.get_ylim()
            xdata, ydata = event.xdata, event.ydata
            if xdata is None or ydata is None: return
            scale = 1.2 if event.step>0 else 0.8
            nw = (cur_xlim[1]-cur_xlim[0])*scale; nh = (cur_ylim[1]-cur_ylim[0])*scale
            rx = (cur_xlim[1]-xdata)/(cur_xlim[1]-cur_xlim[0]); ry = (cur_ylim[1]-ydata)/(cur_ylim[1]-cur_ylim[0])
            ax.set_xlim([xdata-nw*(1-rx), xdata+nw*rx]); ax.set_ylim([ydata-nh*(1-ry), ydata+nh*ry])
            fig.canvas.draw_idle()
        def on_press(event):
            if event.inaxes != ax: return
            if event.button == 1:
                self._rot_active = True
                self._rot_start_elev = ax.elev
                self._rot_start_azim = ax.azim
                self._rot_start_x = event.x
                self._rot_start_y = event.y
            elif event.button == 3:
                self._pan_active = True
                self._pan_start_x = event.xdata
                self._pan_start_y = event.ydata
        def on_release(event):
            self._pan_active = False
            self._rot_active = False
        def on_motion(event):
            if event.inaxes != ax: return
            if getattr(self, '_rot_active', False):
                dx = event.x - self._rot_start_x
                dy = event.y - self._rot_start_y
                ax.view_init(elev=self._rot_start_elev + dy * 0.5,
                             azim=self._rot_start_azim - dx * 0.5)
                fig.canvas.draw_idle()
            elif self._pan_active and event.xdata is not None and event.ydata is not None:
                dx = self._pan_start_x - event.xdata
                dy = self._pan_start_y - event.ydata
                ax.set_xlim(ax.get_xlim() + dx)
                ax.set_ylim(ax.get_ylim() + dy)
                fig.canvas.draw_idle()
        fig.canvas.mpl_connect('scroll_event', zoom_fun)
        fig.canvas.mpl_connect('button_press_event', on_press)
        fig.canvas.mpl_connect('button_release_event', on_release)
        fig.canvas.mpl_connect('motion_notify_event', on_motion)

    def on_tab_changed(self, event):
        if self.notebook.index(self.notebook.select()) == 4:
            self._plot_stpc_graph()

    def show_tree_menu(self, event): self.tree_menu.post(event.x_root, event.y_root)
    def copy_selected_rows(self):
        sel = self.tree.selection()
        if sel:
            rows = ['\t'.join(str(v) for v in self.tree.item(i,'values')) for i in sel]
            self.root.clipboard_clear(); self.root.clipboard_append('\n'.join(rows))
            self.set_status("Copied to clipboard")

if __name__ == "__main__":
    root = tk.Tk(); AllyMindApp(root); root.mainloop()
