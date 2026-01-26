# -*- coding: utf-8 -*-
"""
Created on Sun Nov  2 19:03:22 2025

@author: Emily Gemmill
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path
from aicspylibczi import CziFile
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from tkinter import ttk, messagebox
from scipy.ndimage import rotate

class EggChamberOrientationTool:
    def __init__(self, root_folder_path):
        self.root_folder_path = Path(root_folder_path)
        
        # Find all genotype subfolders with CZI files
        self.genotype_folders = []
        for subfolder in sorted(self.root_folder_path.iterdir()):
            if subfolder.is_dir():
                czi_files = list(subfolder.glob('P*_20x.czi'))
                if czi_files:
                    self.genotype_folders.append(subfolder)
        
        if not self.genotype_folders:
            raise ValueError(f"No subfolders with P*_20x.czi files found in {root_folder_path}")
        
        print(f"Found {len(self.genotype_folders)} genotype folders:")
        for folder in self.genotype_folders:
            print(f"  - {folder.name}")
        
        self.current_folder_idx = 0
        self.current_idx = 0
        self.click_points = []
        self.results = []
        
        # Setup GUI
        self.root = tk.Tk()
        self.root.title("Egg Chamber Orientation Tool")
        self.setup_gui()
        
        # Load first folder
        self.load_current_folder()
        self.load_current_image()
        
    def setup_gui(self):
        """Setup the GUI layout"""
        # Top frame for instructions and controls
        top_frame = ttk.Frame(self.root, padding="10")
        top_frame.grid(row=0, column=0, sticky=(tk.W, tk.E))
        
        # Instructions
        instructions = ttk.Label(top_frame, text=(
            "Instructions:\n"
            "1. First click: Center of egg chamber (or border cell cluster)\n"
            "2. Second click: Lead protrusion (direction of migration)\n"
            "Arrow will show from center → lead protrusion"
        ), font=('Arial', 10))
        instructions.grid(row=0, column=0, columnspan=3, pady=5)
        
        # Progress label
        self.progress_label = ttk.Label(top_frame, text="", font=('Arial', 12, 'bold'))
        self.progress_label.grid(row=1, column=0, columnspan=3, pady=5)
        
        # Genotype label
        self.genotype_label = ttk.Label(top_frame, text="", font=('Arial', 11, 'bold'), foreground='blue')
        self.genotype_label.grid(row=2, column=0, columnspan=3, pady=2)
        
        # Position label
        self.position_label = ttk.Label(top_frame, text="", font=('Arial', 11))
        self.position_label.grid(row=3, column=0, columnspan=3, pady=2)
        
        # Button frame
        button_frame = ttk.Frame(top_frame)
        button_frame.grid(row=4, column=0, columnspan=3, pady=10)
        
        ttk.Button(button_frame, text="← Previous", command=self.prev_image).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear Clicks", command=self.clear_clicks).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Skip", command=self.skip_image).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Next →", command=self.next_image).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Save & Quit", command=self.save_and_quit).pack(side=tk.LEFT, padx=5)
        
        # Canvas frame
        canvas_frame = ttk.Frame(self.root)
        canvas_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create matplotlib figure
        self.fig, (self.ax1, self.ax2) = plt.subplots(1, 2, figsize=(16, 8))
        self.canvas = FigureCanvasTkAgg(self.fig, master=canvas_frame)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Connect click event
        self.canvas.mpl_connect('button_press_event', self.on_click)
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        
    def load_current_folder(self):
        """Load current genotype folder"""
        self.folder_path = self.genotype_folders[self.current_folder_idx]
        self.genotype_name = self.folder_path.name
        self.czi_files = sorted(self.folder_path.glob('P*_20x.czi'))
        self.results = []  # Reset results for new folder
        self.current_idx = 0
        print(f"\n{'='*60}")
        print(f"Loading folder: {self.genotype_name} ({len(self.czi_files)} images)")
        print(f"{'='*60}")
    
    def load_current_image(self):
        """Load and display current image"""
        # Check if we've finished this folder
        if self.current_idx >= len(self.czi_files):
            self.finish_current_folder()
            return
        
        czi_path = self.czi_files[self.current_idx]
        self.position_name = czi_path.stem.replace('_20x', '')
        
        # Update labels
        total_folders = len(self.genotype_folders)
        self.progress_label.config(
            text=f"Folder {self.current_folder_idx + 1}/{total_folders} | Image {self.current_idx + 1}/{len(self.czi_files)}"
        )
        self.genotype_label.config(text=f"Genotype: {self.genotype_name}")
        self.position_label.config(text=f"Position: {self.position_name}")
        
        # Load 20x CZI
        try:
            czi = CziFile(czi_path)
            img_data, _ = czi.read_image()
            self.img_20x = np.squeeze(img_data)
            
            # Create RGB composite
            if self.img_20x.ndim == 3 and self.img_20x.shape[0] <= 4:
                rgb = np.zeros((self.img_20x.shape[1], self.img_20x.shape[2], 3), dtype=np.uint8)
                for c in range(min(3, self.img_20x.shape[0])):
                    channel = self.img_20x[c].astype(np.float32)
                    if channel.max() > channel.min():
                        channel_norm = ((channel - channel.min()) / (channel.max() - channel.min()) * 255)
                        rgb[:, :, c] = channel_norm.astype(np.uint8)
                self.rgb_20x = rgb
            else:
                # Grayscale
                img = self.img_20x.astype(np.float32)
                img = ((img - img.min()) / (img.max() - img.min() + 1e-8) * 255).astype(np.uint8)
                self.rgb_20x = np.stack([img]*3, axis=-1)
        except Exception as e:
            print(f"Error loading 20x image: {e}")
            self.rgb_20x = np.zeros((100, 100, 3), dtype=np.uint8)
        
        # Load 40x mask if exists
        mask_path = self.folder_path / f"{self.position_name}_mask.tif"
        if mask_path.exists():
            try:
                self.mask_40x = np.array(Image.open(mask_path))
            except Exception as e:
                print(f"Error loading mask: {e}")
                self.mask_40x = None
        else:
            self.mask_40x = None
        
        # Reset clicks
        self.click_points = []
        
        # Display
        self.update_display()
        
    def update_display(self):
        """Update the display with current image and clicks"""
        self.ax1.clear()
        self.ax2.clear()
        
        # Display 20x image
        self.ax1.imshow(self.rgb_20x)
        self.ax1.set_title(f'{self.position_name} - 20x Image\nClick: 1) Center, 2) Lead protrusion', 
                          fontsize=12)
        self.ax1.axis('off')
        
        # Draw clicks
        if len(self.click_points) >= 1:
            # Draw centroid
            self.ax1.plot(self.click_points[0][1], self.click_points[0][0], 'g+', 
                         markersize=30, markeredgewidth=4, label='Center')
        
        if len(self.click_points) >= 2:
            # Draw lead point
            self.ax1.plot(self.click_points[1][1], self.click_points[1][0], 'r*', 
                         markersize=30, markeredgewidth=2, label='Lead protrusion')
            
            # Draw arrow
            dy = self.click_points[1][0] - self.click_points[0][0]
            dx = self.click_points[1][1] - self.click_points[0][1]
            self.ax1.arrow(self.click_points[0][1], self.click_points[0][0], 
                          dx, dy, head_width=30, head_length=40, 
                          fc='yellow', ec='yellow', linewidth=4)
            
            # Calculate angle
            angle = np.degrees(np.arctan2(dy, dx))
            self.ax1.text(0.5, 0.95, f'Angle: {angle:.1f}°', 
                         transform=self.ax1.transAxes, fontsize=14, 
                         ha='center', va='top',
                         bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.8))
        
        self.ax1.legend(loc='upper left', fontsize=10)
        
        # Display 40x mask if available
        if self.mask_40x is not None:
            self.ax2.imshow(self.mask_40x, cmap='gray')
            self.ax2.set_title(f'{self.position_name} - 40x Mask', fontsize=12)
        else:
            self.ax2.text(0.5, 0.5, 'No 40x mask available', 
                         transform=self.ax2.transAxes, ha='center', va='center',
                         fontsize=14)
        self.ax2.axis('off')
        
        self.canvas.draw()
        
    def on_click(self, event):
        """Handle mouse clicks"""
        if event.inaxes != self.ax1:
            return
        
        if len(self.click_points) >= 2:
            messagebox.showinfo("Info", "Already have 2 points. Click 'Clear Clicks' to start over or 'Next' to continue.")
            return
        
        # Add click point (y, x)
        self.click_points.append([event.ydata, event.xdata])
        
        if len(self.click_points) == 2:
            # Calculate angle
            dy = self.click_points[1][0] - self.click_points[0][0]
            dx = self.click_points[1][1] - self.click_points[0][1]
            angle = np.degrees(np.arctan2(dy, dx))
            print(f"{self.position_name}: Angle = {angle:.1f}°")
        
        self.update_display()
        
    def clear_clicks(self):
        """Clear current clicks"""
        self.click_points = []
        self.update_display()
        
    def skip_image(self):
        """Skip current image (save as NaN)"""
        self.results.append({
            'genotype': self.genotype_name,
            'position': self.position_name,
            'angle': np.nan,
            'skipped': True
        })
        self.current_idx += 1
        
        # Check if we've reached the end
        if self.current_idx >= len(self.czi_files):
            self.finish_current_folder()
        else:
            self.load_current_image()
        
    def next_image(self):
        """Save current and move to next"""
        if len(self.click_points) < 2:
            response = messagebox.askyesno("Incomplete", 
                                          "You haven't marked both points. Skip this image?")
            if response:
                self.skip_image()
            return
        
        # Calculate angle
        dy = self.click_points[1][0] - self.click_points[0][0]
        dx = self.click_points[1][1] - self.click_points[0][1]
        angle = np.degrees(np.arctan2(dy, dx))
        
        # Save result
        self.results.append({
            'genotype': self.genotype_name,
            'position': self.position_name,
            'angle': angle,
            'center_y': self.click_points[0][0],
            'center_x': self.click_points[0][1],
            'lead_y': self.click_points[1][0],
            'lead_x': self.click_points[1][1],
            'skipped': False
        })
        
        self.current_idx += 1
        
        # Check if we've reached the end of this folder
        if self.current_idx >= len(self.czi_files):
            self.finish_current_folder()
        else:
            self.load_current_image()
    
    def finish_current_folder(self):
        """Save results for current folder and move to next folder or quit"""
        # Save results for this folder
        if self.results:
            df = pd.DataFrame(self.results)
            output_path = self.folder_path / 'orientation_angles.csv'
            df.to_csv(output_path, index=False)
            print(f"Saved {len(self.results)} annotations to: {output_path}")
        
        # Check if there are more folders
        if self.current_folder_idx + 1 < len(self.genotype_folders):
            response = messagebox.askyesno(
                "Folder Complete", 
                f"Finished {self.genotype_name}!\n\n"
                f"Saved {len(self.results)} annotations.\n\n"
                f"Continue to next genotype folder?"
            )
            
            if response:
                self.current_folder_idx += 1
                self.load_current_folder()
                self.load_current_image()
            else:
                self.save_and_quit()
        else:
            # All folders complete
            messagebox.showinfo(
                "All Complete!", 
                f"Finished all {len(self.genotype_folders)} genotype folders!\n\n"
                f"Results saved in each folder's orientation_angles.csv"
            )
            self.root.quit()
        
    def prev_image(self):
        """Go back to previous image"""
        if self.current_idx > 0:
            self.current_idx -= 1
            # Remove last result if it exists
            if self.results and self.results[-1]['position'] == self.czi_files[self.current_idx].stem.replace('_20x', ''):
                self.results.pop()
            self.load_current_image()
        
    def save_and_quit(self):
        """Save current results and quit"""
        if self.results:
            df = pd.DataFrame(self.results)
            output_path = self.folder_path / 'orientation_angles.csv'
            df.to_csv(output_path, index=False)
            print(f"Saved {len(self.results)} annotations to: {output_path}")
            messagebox.showinfo("Saved", 
                              f"Saved {len(self.results)} annotations for {self.genotype_name}")
        
        self.root.quit()
        
    def run(self):
        """Run the GUI"""
        self.root.mainloop()

if __name__ == "__main__":
    root_folder = input("Enter root folder path (containing genotype subfolders): ").strip().strip('"').strip("'")
    
    try:
        tool = EggChamberOrientationTool(root_folder)
        tool.run()
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()