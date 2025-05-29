"""
# ==============================================================================
# Author:       Victor S Caricatte De Araújo
# Email:        victorleniwys@gmail.com or victorsc@ufmg.br
# Intitution:   Universidade federal de Minas Gerais
# Version:      1.1
# Date:         Abr, 14
# ...................................
# ==============================================================================
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import os
from Bio.Seq import Seq
import webbrowser
from tkinter.scrolledtext import ScrolledText

class DNAPalindromeFinder:
    """
    Application for identifying palindromic sequences in DNA.
    
    DNA palindromes are sequences that are identical to their reverse complement.
    These sequences are important as they often represent binding sites for
    restriction enzymes.
    """
    
    def __init__(self, root):
        """Initialize the application with the main window"""
        self.root = root
        self._setup_main_window()
        self.create_widgets()
        self._set_default_values()
    
    def _setup_main_window(self):
        """Configure main window properties"""
        self.root.title("DNA Palindrome Finder - Bioinformatics Tool")
        self.root.geometry("900x700")
        self.root.minsize(800, 600)
        self.root.option_add('*tearOff', False)
        
        # Configure style
        style = ttk.Style()
        style.configure('TButton', padding=5)
        style.configure('TEntry', padding=5)
        style.configure('TNotebook.Tab', padding=[10, 5])
    
    def _set_default_values(self):
        """Set default values for fields"""
        self.min_length.set("4")
        self.max_length.set("12")
        self.status_var.set("Ready to analyze DNA sequences")
    
    def create_widgets(self):
        """Create all interface elements"""
        # Create notebook (tabs)
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Main tab (analysis)
        self.main_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.main_frame, text="DNA Analysis")
        
        # Help tab
        self.help_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.help_frame, text="Help")
        
        # Configure main tab
        self._setup_main_tab()
        
        # Configure help tab
        self._setup_help_tab()
        
        # Create menu
        self._create_menu()
        
        # Status bar
        self.status_var = tk.StringVar()
        self.status_bar = ttk.Label(self.root, textvariable=self.status_var, relief=tk.SUNKEN)
        self.status_bar.pack(side=tk.BOTTOM, fill=tk.X)
    
    def _setup_main_tab(self):
        """Configure the main analysis tab"""
        # Input frame
        input_frame = ttk.LabelFrame(self.main_frame, text="DNA Input", padding="10")
        input_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Text area for DNA sequence
        self.dna_entry = ScrolledText(input_frame, height=10, wrap=tk.WORD, font=('Courier', 10))
        self.dna_entry.pack(fill=tk.BOTH, expand=True)
        
        # Settings frame
        settings_frame = ttk.Frame(self.main_frame)
        settings_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Minimum length
        ttk.Label(settings_frame, text="Min length:").pack(side=tk.LEFT)
        self.min_length = tk.StringVar()
        ttk.Entry(settings_frame, textvariable=self.min_length, width=5).pack(side=tk.LEFT, padx=(0, 20))
        
        # Maximum length
        ttk.Label(settings_frame, text="Max length:").pack(side=tk.LEFT)
        self.max_length = tk.StringVar()
        ttk.Entry(settings_frame, textvariable=self.max_length, width=5).pack(side=tk.LEFT)
        
        # Analyze button
        ttk.Button(settings_frame, text="Analyze DNA", command=self.find_palindromes).pack(side=tk.RIGHT)
        
        # Results frame
        results_frame = ttk.LabelFrame(self.main_frame, text="Results", padding="10")
        results_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Treeview for results
        self.results_tree = ttk.Treeview(results_frame, columns=('Position', 'Sequence', 'Length'), show='headings')
        self.results_tree.heading('Position', text='Position')
        self.results_tree.heading('Sequence', text='Sequence')
        self.results_tree.heading('Length', text='Length')
        self.results_tree.column('Position', width=100, anchor=tk.CENTER)
        self.results_tree.column('Sequence', width=200, anchor=tk.CENTER)
        self.results_tree.column('Length', width=100, anchor=tk.CENTER)
        
        # Scrollbar and packing
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.results_tree.yview)
        self.results_tree.configure(yscroll=scrollbar.set)
        
        # Action buttons
        button_frame = ttk.Frame(results_frame)
        button_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=(5, 0))
        
        ttk.Button(button_frame, text="Copy Results", command=self.copy_results).pack(side=tk.LEFT)
        ttk.Button(button_frame, text="Save Results", command=self.save_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear All", command=self.clear_all).pack(side=tk.RIGHT)
        
        # Pack treeview and scrollbar
        self.results_tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    
    def _setup_help_tab(self):
        """Configure the help/documentation tab"""
        # Main frame
        help_content = ScrolledText(self.help_frame, wrap=tk.WORD, padx=10, pady=10)
        help_content.pack(fill=tk.BOTH, expand=True)
        
        # Help content
        help_text = """=== DNA Palindrome Finder ===

DESCRIPTION:
This program identifies palindromic sequences in DNA, which are sequences that are identical when compared to their reverse complement strand. These sequences are important in molecular biology as they often represent binding sites for restriction enzymes.

HOW TO USE:
1. Enter your DNA sequence in the text area (only A, T, C, G characters)
2. Set the minimum and maximum length of palindromes to search for
3. Click "Analyze DNA" to start the search
4. Results will be displayed in the table below

FEATURES:
- Load sequences from FASTA files
- Save results in text or CSV format
- Copy results to clipboard
- Configure palindrome length parameters

TECHNICAL DEFINITION:
A DNA palindrome is a nucleotide sequence that is identical to its reverse complement strand. For example:
Sequence: 5'-GAATTC-3'
Reverse complement: 3'-CTTAAG-5'
When reversed: 5'-GAATTC-3'

This is the recognition site for the EcoRI restriction enzyme.

MENU:
- File: Load/save files and exit the program
- Help: Documentation and program information

DEVELOPED BY:
Victor Silveira Caricatte/ Universidade Federal de Minas Gerais.

VERSION: 1.1 (Python 3.x)
"""
        help_content.insert(tk.END, help_text)
        help_content.config(state=tk.DISABLED)
        
    
    def _create_menu(self):
        """Create the menu bar"""
        menubar = tk.Menu(self.root)
        self.root.config(menu=menubar)
        
        # File Menu
        file_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="Open FASTA file...", command=self.load_fasta_file)
        file_menu.add_command(label="Save results...", command=self.save_results)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        
        # Help Menu
        help_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Help", menu=help_menu)
        help_menu.add_command(label="Documentation", command=self.show_documentation)
        help_menu.add_command(label="About", command=self.show_about)
    
    def validate_dna_sequence(self, sequence):
        """Validate that the sequence contains only valid bases (A, T, C, G)"""
        sequence = sequence.upper().replace("\n", "").replace(" ", "")
        if not sequence:
            raise ValueError("DNA sequence is empty")
        invalid_bases = set(sequence) - {'A', 'T', 'C', 'G'}
        if invalid_bases:
            raise ValueError(f"Invalid bases found: {', '.join(invalid_bases)}")
        return sequence
    
    def is_palindrome(self, substring):
        """Efficiently checks if a sequence is a palindrome"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return substring == ''.join([complement[base] for base in substring[::-1]])
    
    def find_dna_palindromes(self, sequence, min_len, max_len):
        """Optimized version of palindrome search"""
        palindromes = []
        n = len(sequence)
        
        for length in range(min_len, min(max_len, n) + 1):
            for i in range(n - length + 1):
                substring = sequence[i:i+length]
                if self.is_palindrome(substring):
                    palindromes.append((i+1, substring))  # +1 for 1-based position
        
        return sorted(palindromes, key=lambda x: x[0])  # Sort by position
    
    def find_palindromes(self):
        """Controls the palindrome search process"""
        try:
            # Clear previous results
            self.results_tree.delete(*self.results_tree.get_children())
            
            # Get and validate sequence
            dna_sequence = self.validate_dna_sequence(self.dna_entry.get("1.0", tk.END))
            
            # Get parameters
            min_len = max(2, int(self.min_length.get()))  # Ensure minimum of 2
            max_len = max(min_len, int(self.max_length.get()))  # Ensure max >= min
            
            # Update interface during long processing
            self.status_var.set("Searching for palindromes...")
            self.root.update()
            
            # Search palindromes in a separate thread to avoid freezing the interface
            import threading
            threading.Thread(
                target=self._find_palindromes_thread,
                args=(dna_sequence, min_len, max_len),
                daemon=True
            ).start()
            
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            self.status_var.set("Input data error")
        except Exception as e:
            messagebox.showerror("Error", f"Unexpected error: {str(e)}")
            self.status_var.set("Processing error")
    
    def _find_palindromes_thread(self, sequence, min_len, max_len):
        """Runs in separate thread to avoid freezing the GUI"""
        try:
            palindromes = self.find_dna_palindromes(sequence, min_len, max_len)
            
            # Update interface in main thread
            self.root.after(0, self._display_results, palindromes)
            
        except Exception as e:
            self.root.after(0, messagebox.showerror, "Error", f"Error during search:\n{str(e)}")
            self.root.after(0, lambda: self.status_var.set("Processing error"))
    
    def _display_results(self, palindromes):
        """Display results in the interface"""
        if not palindromes:
            messagebox.showinfo("Information", "No palindromes found with the specified criteria.")
            self.status_var.set("Search completed - No palindromes found")
            return
        
        for pos, seq in palindromes:
            self.results_tree.insert('', tk.END, values=(f"{pos}-{pos+len(seq)-1}", seq, len(seq)))
        
        self.status_var.set(f"Search completed - {len(palindromes)} palindromes found")
    
    def load_fasta_file(self):
        """Load sequence from a FASTA file"""
        filepath = filedialog.askopenfilename(
            filetypes=[("FASTA files", "*.fasta;*.fa;*.faa;*.fna"), ("Text files", "*.txt"), ("All files", "*.*")],
            title="Select a FASTA file"
        )
        if filepath:
            try:
                with open(filepath, 'r') as file:
                    content = file.read()
                    # Extract sequence (ignore header and comment lines)
                    sequence = ''.join(line.strip() for line in content.split('\n') if not line.startswith('>'))
                    self.dna_entry.delete("1.0", tk.END)
                    self.dna_entry.insert("1.0", sequence)
                    self.status_var.set(f"File loaded: {os.path.basename(filepath)}")
            except Exception as e:
                messagebox.showerror("Error", f"Could not read file:\n{str(e)}")
                self.status_var.set("Error loading file")
    
    def save_results(self):
        """Save results to a file"""
        if not self.results_tree.get_children():
            messagebox.showwarning("Warning", "No results to save")
            return
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")],
            title="Save results as"
        )
        
        if filepath:
            try:
                with open(filepath, 'w') as file:
                    # Write header
                    file.write("Position\tSequence\tLength\n")
                    # Write data
                    for item in self.results_tree.get_children():
                        values = self.results_tree.item(item, 'values')
                        file.write("\t".join(values) + "\n")
                self.status_var.set(f"Results saved to {os.path.basename(filepath)}")
            except Exception as e:
                messagebox.showerror("Error", f"Could not save file:\n{str(e)}")
                self.status_var.set("Error saving results")
    
    def copy_results(self):
        """Copy results to clipboard"""
        if not self.results_tree.get_children():
            messagebox.showwarning("Warning", "No results to copy")
            return
        
        result_text = "Position\tSequence\tLength\n"
        for item in self.results_tree.get_children():
            values = self.results_tree.item(item, 'values')
            result_text += "\t".join(values) + "\n"
        
        self.root.clipboard_clear()
        self.root.clipboard_append(result_text)
        self.status_var.set("Results copied to clipboard")
    
    def clear_all(self):
        """Clear all fields"""
        self.dna_entry.delete("1.0", tk.END)
        self.results_tree.delete(*self.results_tree.get_children())
        self.min_length.set("4")
        self.max_length.set("12")
        self.status_var.set("Fields cleared - Ready for new analysis")
    
    def show_documentation(self):
        """Show the help tab"""
        self.notebook.select(self.help_frame)
    
    def show_about(self):
        """Show the 'About' window"""
        about_window = tk.Toplevel(self.root)
        about_window.title("About DNA Palindrome Finder")
        about_window.geometry("400x300")
        
        ttk.Label(about_window, 
                 text="DNA Palindrome Finder\n\nVersion 1.1\n\n"
                      "Bioinformatics Tool\n\n"
                      "Developed for DNA sequence analysis\n"
                      "and identification of palindromic\n"
                      "restriction sites",
                 justify=tk.CENTER).pack(pady=20)
        
        ttk.Button(about_window, text="Close", command=about_window.destroy).pack(pady=10)

if __name__ == "__main__":
    root = tk.Tk()
    app = DNAPalindromeFinder(root)
    root.mainloop()
