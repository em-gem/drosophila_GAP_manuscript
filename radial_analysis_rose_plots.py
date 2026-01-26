# -*- coding: utf-8 -*-
"""
Created on Fri Nov  7 14:56:33 2025

@author: Emily Gemmill
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from skimage import io, measure
from scipy.signal import find_peaks
import pandas as pd

def extract_radial_profile(binary_img, n_bins=360, pad=5):
    """    Extract radial profile from binary mask    """
    binary_img = np.pad(binary_img, pad_width=pad, mode='constant', constant_values=0)

    coords = np.argwhere(binary_img > 0)
    if coords.size == 0:
        return np.array([]), np.array([])

    cy, cx = coords.mean(axis=0)

    contours = measure.find_contours(binary_img, 0.5)
    if len(contours) == 0:
        return np.array([]), np.array([])

    contour = max(contours, key=len)
    y, x = contour[:, 0] - cy, contour[:, 1] - cx
    angles = np.arctan2(y, x) % (2*np.pi)
    radii = np.sqrt(x**2 + y**2)

    theta_grid = np.linspace(0, 2*np.pi, n_bins, endpoint=False)
    r_grid = np.interp(theta_grid, np.sort(angles), radii[np.argsort(angles)], period=2*np.pi)
    return theta_grid, r_grid

def find_protrusions(theta_grid, r_grid, mean_radius):
    """Find protrusions in radial profile"""
    if len(r_grid) == 0:
        return np.array([], dtype=int)
    height_thresh = mean_radius * 1.20
    prominence_thresh = mean_radius * 0.5
    peaks, _ = find_peaks(r_grid, height=height_thresh, prominence=prominence_thresh)
    return peaks

def filter_ectopic_protrusions(protrusion_angles, exclusion_zone_deg=20):
    """
    Filter out lead protrusions (those within ±exclusion_zone_deg of 0°).
    
    Args:
        protrusion_angles: array of angles in radians (oriented, so lead is at 0°)
        exclusion_zone_deg: degrees to exclude around 0° (default 20)
    
    Returns:
        array of ectopic (non-lead) protrusion angles
    """
    if len(protrusion_angles) == 0:
        return np.array([])
    
    exclusion_rad = np.deg2rad(exclusion_zone_deg)
    
    # Keep angles that are NOT within ±exclusion_zone of 0° (or 2π)
    # Account for wraparound at 0/2π
    ectopic_mask = (protrusion_angles > exclusion_rad) & (protrusion_angles < (2*np.pi - exclusion_rad))
    
    return protrusion_angles[ectopic_mask]

def apply_orientation(angles, orientation_angle_deg):
    """
    Apply orientation correction to angles.
    Subtracts the orientation angle to align all images with oocyte at 0°.
    
    Args:
        angles: numpy array of angles in radians
        orientation_angle_deg: orientation angle in degrees (from orientation_angles.csv)
    
    Returns:
        rotated angles in radians, normalized to [0, 2π)
    """
    orientation_rad = np.deg2rad(orientation_angle_deg)
    rotated = (angles - orientation_rad) % (2 * np.pi)
    return rotated

def plot_rose_triple(theta_orig, r_grid_orig, peaks_orig, 
                     theta_oriented, r_grid_oriented, peaks_oriented,
                     peaks_ectopic, save_path, title="Rose Plot Comparison",
                     exclusion_zone_deg=20):
    """Plot three-panel comparison: original, oriented, and ectopic only"""
    if len(theta_orig) == 0 or len(r_grid_orig) == 0:
        print(f"Skipping plot {save_path} (empty data)")
        return
    
    fig = plt.figure(figsize=(18, 5))
    
    # Original (unoriented)
    ax1 = plt.subplot(1, 3, 1, projection='polar')
    ax1.plot(theta_orig, r_grid_orig, 'b-', lw=1.5)
    if len(peaks_orig) > 0:
        ax1.scatter(theta_orig[peaks_orig], r_grid_orig[peaks_orig], 
                   c="r", s=50, zorder=5, label="Protrusions")
        ax1.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=9)
    ax1.set_title("Original (Unoriented)", pad=20)
    ax1.set_theta_zero_location('E')
    ax1.set_rlabel_position(135)  # Move radial labels away from data
    
    # Oriented (all protrusions)
    ax2 = plt.subplot(1, 3, 2, projection='polar')
    ax2.plot(theta_oriented, r_grid_oriented, 'b-', lw=1.5)
    if len(peaks_oriented) > 0:
        ax2.scatter(theta_oriented[peaks_oriented], r_grid_oriented[peaks_oriented], 
                   c="r", s=50, zorder=5, label="All protrusions")
        ax2.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=9)
    ax2.plot([0, 0], [0, max(r_grid_oriented)], 'g--', lw=2, alpha=0.5)
    ax2.set_title("Oriented (Lead at 0°)", pad=20)
    ax2.set_theta_zero_location('E')
    ax2.set_rlabel_position(135)
    
    # Ectopic only (exclusion zone marked)
    ax3 = plt.subplot(1, 3, 3, projection='polar')
    ax3.plot(theta_oriented, r_grid_oriented, 'b-', lw=1.5, alpha=0.3)
    if len(peaks_ectopic) > 0:
        ax3.scatter(theta_oriented[peaks_ectopic], r_grid_oriented[peaks_ectopic], 
                   c="orange", s=50, zorder=5, label="Ectopic protrusions")
        ax3.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=9)
    
    # Mark exclusion zone
    exclusion_rad = np.deg2rad(exclusion_zone_deg)
    ax3.fill_between([0, exclusion_rad], 0, max(r_grid_oriented), 
                     alpha=0.2, color='gray', label=f'Excluded (±{exclusion_zone_deg}°)')
    ax3.fill_between([2*np.pi - exclusion_rad, 2*np.pi], 0, max(r_grid_oriented), 
                     alpha=0.2, color='gray')
    ax3.plot([0, 0], [0, max(r_grid_oriented)], 'r--', lw=2, alpha=0.5)
    ax3.set_title(f"Ectopic Only (Lead ±{exclusion_zone_deg}° excluded)", pad=20)
    ax3.set_theta_zero_location('E')
    ax3.set_rlabel_position(135)
    
    plt.suptitle(title, y=1.02, fontsize=12)
    plt.tight_layout()
    plt.savefig(save_path, dpi=600, bbox_inches='tight')
    plt.close(fig)

def rayleigh_test(angles):
    """Perform Rayleigh test for circular uniformity"""
    angles = np.array(angles, dtype=float)
    n = len(angles)
    if n == 0:
        return np.nan, np.nan
    R = np.abs(np.sum(np.exp(1j * angles))) / n
    z = n * R**2
    p = np.exp(-z) * (1 + (2*z - z**2) / (4*n) -
                      (24*z - 132*z**2 + 76*z**3 - 9*z**4) / (288*n**2))
    return R, p

def process_mask_with_orientation(genotype, fname, img_path, orientation_angle, 
                                   threshold_mode, global_mean_radius, geno_out, 
                                   exclusion_zone_deg=20):
    """Process mask with orientation correction"""
    mask = io.imread(img_path) > 0
    theta_grid, r_grid = extract_radial_profile(mask)
    if len(r_grid) == 0:
        print(f"Skipping {fname} (empty mask or contour failed)")
        return None

    mean_radius = np.mean(r_grid)

    # Find protrusions
    if threshold_mode == "individual":
        peaks = find_protrusions(theta_grid, r_grid, mean_radius)
    else:
        peaks = find_protrusions(theta_grid, r_grid, global_mean_radius)

    # Get original (unoriented) protrusion angles
    protrusion_angles_orig = theta_grid[peaks] if len(peaks) > 0 else np.array([])
    
    # Calculate Rayleigh statistics for individual cluster
    R_orig, p_orig = rayleigh_test(protrusion_angles_orig)
    
    # Apply orientation if available
    if not np.isnan(orientation_angle):
        # Rotate theta grid to orient egg chamber
        theta_oriented = apply_orientation(theta_grid, orientation_angle)
        protrusion_angles_oriented = theta_oriented[peaks] if len(peaks) > 0 else np.array([])
        
        # Calculate Rayleigh for oriented angles
        R_oriented, p_oriented = rayleigh_test(protrusion_angles_oriented)
        
        # Filter ectopic protrusions (exclude lead ±20°)
        ectopic_protrusion_angles = filter_ectopic_protrusions(protrusion_angles_oriented, exclusion_zone_deg)
        
        # Calculate Rayleigh for ectopic angles
        R_ectopic, p_ectopic = rayleigh_test(ectopic_protrusion_angles)
        
        # Find which peaks are ectopic
        if len(peaks) > 0 and len(ectopic_protrusion_angles) > 0:
            # Map ectopic angles back to peak indices
            ectopic_peaks = []
            for angle in ectopic_protrusion_angles:
                # Find the peak index that corresponds to this angle
                idx = np.argmin(np.abs(protrusion_angles_oriented - angle))
                if np.abs(protrusion_angles_oriented[idx] - angle) < 1e-6:  # Match found
                    ectopic_peaks.append(peaks[idx])
            ectopic_peaks = np.array(ectopic_peaks)
        else:
            ectopic_peaks = np.array([])
        
        # Save three-panel comparison plot
        save_path = os.path.join(geno_out, fname.replace(".tif", "_rose_comparison.png"))
        plot_rose_triple(theta_grid, r_grid, peaks,
                        theta_oriented, r_grid, peaks,
                        ectopic_peaks, save_path, 
                        title=f"{genotype} - {fname}",
                        exclusion_zone_deg=exclusion_zone_deg)
    else:
        print(f"  No orientation data for {fname}, skipping oriented plot")
        theta_oriented = theta_grid
        protrusion_angles_oriented = protrusion_angles_orig
        ectopic_protrusion_angles = np.array([])
        R_oriented, p_oriented = np.nan, np.nan
        R_ectopic, p_ectopic = np.nan, np.nan

    return {
        "genotype": genotype,
        "file": fname,
        "n_protrusions": len(peaks),
        "n_ectopic_protrusions": len(ectopic_protrusion_angles),
        "mean_radius": mean_radius,
        "orientation_angle": orientation_angle,
        "protrusion_angles_original": ";".join(map(str, protrusion_angles_orig)),
        "protrusion_angles_oriented": ";".join(map(str, protrusion_angles_oriented)),
        "protrusion_angles_ectopic": ";".join(map(str, ectopic_protrusion_angles)),
        "protrusion_radii": ";".join(map(str, r_grid[peaks])),
        "rayleigh_R_original": R_orig,
        "rayleigh_p_original": p_orig,
        "rayleigh_R_oriented": R_oriented,
        "rayleigh_p_oriented": p_oriented,
        "rayleigh_R_ectopic": R_ectopic,
        "rayleigh_p_ectopic": p_ectopic
    }

def load_orientation_data(genotype_path):
    """Load orientation angles from CSV if available"""
    csv_path = os.path.join(genotype_path, 'orientation_angles.csv')
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        # Create dictionary mapping position to angle
        orientation_dict = {}
        for _, row in df.iterrows():
            if not row['skipped']:
                orientation_dict[row['position']] = row['angle']
        return orientation_dict
    return {}

def analyze_dataset_with_orientation(root_dir, output_dir, mode="individual"):
    """Analyze dataset with orientation correction"""
    results = []
    os.makedirs(output_dir, exist_ok=True)

    # Calculate global mean radius if needed
    global_mean_radius = None
    if mode == "global":
        global_radii = []
        for genotype in os.listdir(root_dir):
            geno_path = os.path.join(root_dir, genotype)
            if not os.path.isdir(geno_path): 
                continue
            for fname in os.listdir(geno_path):
                if fname.endswith("_mask.tif"):
                    img = io.imread(os.path.join(geno_path, fname)) > 0
                    _, r_grid = extract_radial_profile(img)
                    if len(r_grid) > 0:
                        global_radii.append(np.mean(r_grid))
        if len(global_radii) == 0:
            raise ValueError("No valid masks found to compute global mean radius.")
        global_mean_radius = np.mean(global_radii)
        print(f"Global mean radius: {global_mean_radius:.2f} pixels")

    summaries = []
    for genotype in os.listdir(root_dir):
        geno_path = os.path.join(root_dir, genotype)
        if not os.path.isdir(geno_path): 
            continue
        geno_out = os.path.join(output_dir, genotype)
        os.makedirs(geno_out, exist_ok=True)

        print(f"\nProcessing genotype: {genotype}")
        
        # Load orientation data
        orientation_dict = load_orientation_data(geno_path)
        print(f"  Found orientation data for {len(orientation_dict)} positions")

        genotype_results = []
        for fname in os.listdir(geno_path):
            if fname.endswith("_mask.tif"):
                position = fname.replace("_mask.tif", "")
                orientation_angle = orientation_dict.get(position, np.nan)
                
                img_path = os.path.join(geno_path, fname)
                res = process_mask_with_orientation(genotype, fname, img_path, 
                                                   orientation_angle, mode, 
                                                   global_mean_radius, geno_out)
                if res is not None:
                    genotype_results.append(res)

        if len(genotype_results) == 0:
            print(f"  No valid masks found for genotype {genotype}")
            continue

        results.extend(genotype_results)

        # Collect angles for composite plots
        all_angles_orig = []
        all_angles_oriented = []
        all_angles_ectopic = []
        
        for r in genotype_results:
            if r["protrusion_angles_original"]:
                all_angles_orig.extend(map(float, r["protrusion_angles_original"].split(";")))
            if r["protrusion_angles_oriented"]:
                all_angles_oriented.extend(map(float, r["protrusion_angles_oriented"].split(";")))
            if r["protrusion_angles_ectopic"]:
                all_angles_ectopic.extend(map(float, r["protrusion_angles_ectopic"].split(";")))

        if len(all_angles_orig) == 0:
            print(f"  No protrusions detected for genotype {genotype}")
            continue

        # Make composite plots (normalized by number of clusters)
        n_bins = 36
        n_clusters = len(genotype_results)
        
        # Original composite
        hist_orig, bin_edges = np.histogram(all_angles_orig, bins=n_bins, range=(0, 2*np.pi))
        hist_orig_norm = hist_orig / n_clusters  # Normalize by number of clusters
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        
        # Oriented
        if len(all_angles_oriented) > 0:
            hist_oriented, _ = np.histogram(all_angles_oriented, bins=n_bins, range=(0, 2*np.pi))
            hist_oriented_norm = hist_oriented / n_clusters
        
        # Ectopic
        if len(all_angles_ectopic) > 0:
            hist_ectopic, _ = np.histogram(all_angles_ectopic, bins=n_bins, range=(0, 2*np.pi))
            hist_ectopic_norm = hist_ectopic / n_clusters
        
        fig = plt.figure(figsize=(18, 5))
        
        # Original (unoriented)
        ax1 = plt.subplot(1, 3, 1, projection='polar')
        ax1.bar(bin_centers, hist_orig_norm, width=2*np.pi/n_bins, alpha=0.6, align="center", color='steelblue')
        ax1.set_title(f"Original (Unoriented)\n{genotype}\n({n_clusters} clusters)", pad=20)
        ax1.set_ylabel('Protrusions per cluster', labelpad=30)
        ax1.set_theta_zero_location('E')
        ax1.set_rlabel_position(135)
        
        # Oriented (all protrusions)
        ax2 = plt.subplot(1, 3, 2, projection='polar')
        if len(all_angles_oriented) > 0:
            ax2.bar(bin_centers, hist_oriented_norm, width=2*np.pi/n_bins, alpha=0.6, align="center", color='coral')
            max_height = max(hist_oriented_norm)
            ax2.plot([0, 0], [0, max_height], 'g--', lw=3, alpha=0.7, label='Lead (0°)')
            ax2.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=9)
        ax2.set_title(f"Oriented - All Protrusions\n{genotype}\n({n_clusters} clusters)", pad=20)
        ax2.set_ylabel('Protrusions per cluster', labelpad=30)
        ax2.set_theta_zero_location('E')
        ax2.set_rlabel_position(135)
        
        # Ectopic only
        ax3 = plt.subplot(1, 3, 3, projection='polar')
        if len(all_angles_ectopic) > 0:
            ax3.bar(bin_centers, hist_ectopic_norm, width=2*np.pi/n_bins, alpha=0.6, align="center", color='orange')
            
            # Mark exclusion zone
            exclusion_rad = np.deg2rad(20)
            max_height = max(hist_ectopic_norm) if len(hist_ectopic_norm) > 0 else 1
            ax3.fill_between([0, exclusion_rad], 0, max_height, alpha=0.2, color='gray')
            ax3.fill_between([2*np.pi - exclusion_rad, 2*np.pi], 0, max_height, alpha=0.2, color='gray')
            ax3.plot([0, 0], [0, max_height], 'r--', lw=3, alpha=0.7, label='Lead excluded (±20°)')
            ax3.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=9)
        ax3.set_title(f"Ectopic Only (Lead ±20° excluded)\n{genotype}\n({n_clusters} clusters)", pad=20)
        ax3.set_ylabel('Protrusions per cluster', labelpad=30)
        ax3.set_theta_zero_location('E')
        ax3.set_rlabel_position(135)
        
        plt.suptitle(f"Composite Rose Plots (Normalized) - {genotype}", y=1.02, fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(geno_out, f"{genotype}_composite_comparison_normalized.png"), dpi=600, bbox_inches='tight')
        plt.close(fig)

        # Summary stats
        n_imgs = len(genotype_results)
        n_protrusions = [r["n_protrusions"] for r in genotype_results]
        n_ectopic = [r["n_ectopic_protrusions"] for r in genotype_results]
        n_oriented = sum(1 for r in genotype_results if not np.isnan(r["orientation_angle"]))
        
        R_orig, p_orig = rayleigh_test(all_angles_orig)
        R_oriented, p_oriented = rayleigh_test(all_angles_oriented) if len(all_angles_oriented) > 0 else (np.nan, np.nan)
        R_ectopic, p_ectopic = rayleigh_test(all_angles_ectopic) if len(all_angles_ectopic) > 0 else (np.nan, np.nan)

        summaries.append({
            "genotype": genotype,
            "n_images": n_imgs,
            "n_oriented": n_oriented,
            "mean_n_protrusions": np.mean(n_protrusions),
            "mean_n_ectopic_protrusions": np.mean(n_ectopic),
            "prop_zero": np.mean([n==0 for n in n_protrusions]),
            "prop_one": np.mean([n==1 for n in n_protrusions]),
            "prop_two_or_more": np.mean([n>=2 for n in n_protrusions]),
            "mean_radius": np.mean([r["mean_radius"] for r in genotype_results]),
            "rayleigh_R_original": R_orig,
            "rayleigh_p_original": p_orig,
            "rayleigh_R_oriented": R_oriented,
            "rayleigh_p_oriented": p_oriented,
            "rayleigh_R_ectopic": R_ectopic,
            "rayleigh_p_ectopic": p_ectopic
        })

    # Save detailed results (now includes per-cluster Rayleigh statistics)
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(output_dir, f"protrusion_data_with_orientation_{mode}.csv"), index=False)

    # Save summary
    df_summary = pd.DataFrame(summaries)
    df_summary.to_csv(os.path.join(output_dir, f"summary_with_orientation_{mode}.csv"), index=False)

    print(f"\n{'='*60}")
    print("Analysis complete!")
    print(f"Results saved to: {output_dir}")
    print(f"{'='*60}")

    return df, df_summary

# Example usage
if __name__ == "__main__":
    root_dir = input("Enter root directory (containing genotype subfolders): ").strip().strip('"').strip("'")
    output_dir = input("Enter output directory: ").strip().strip('"').strip("'")
    
    print("\nRunning analysis with orientation correction...")
    print("="*60)
    
    df_individual, summary_individual = analyze_dataset_with_orientation(
        root_dir, 
        os.path.join(output_dir, "individual_oriented"), 
        mode="individual"
    )
    
    df_global, summary_global = analyze_dataset_with_orientation(
        root_dir, 
        os.path.join(output_dir, "global_oriented"), 
        mode="global"
    )
    
    print("\nDone! Check the output directory for:")
    print("  - Individual rose plot comparisons (3 panels: original, oriented, ectopic)")
    print("  - Composite rose plot comparisons (3 panels per genotype)")
    print("  - CSV files with original, oriented, and ectopic angles")
    print("  - Per-cluster Rayleigh statistics in detailed CSV")
    print("  - Summary statistics with Rayleigh tests for all three analyses")