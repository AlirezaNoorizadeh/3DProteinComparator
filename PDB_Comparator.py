import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def read_pdb(file_path):
    """خواندن فایل PDB و استخراج مختصات کربن آلفا"""
    ca_coords = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    atom_name = line[12:16].strip()
                    if atom_name == "CA":  # فقط کربن آلفا
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ca_coords.append([x, y, z])
        return np.array(ca_coords)
    except FileNotFoundError:
        print(f"Error: File {file_path} not found!")
        return None
    except Exception as e:
        print(f"Error reading file {file_path}: {str(e)}")
        return None

def kabsch_algorithm(P, Q):
    """الگوریتم Kabsch برای همترازی دو مجموعه نقطه P و Q"""
    # مرکزگرایی
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    # محاسبه ماتریس همبستگی
    H = P_centered.T @ Q_centered

    # تجزیه مقدار تکین (SVD)
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # بررسی reflection (دترمینان باید ۱+ باشد)
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    # اعمال چرخش و انتقال به Q
    Q_aligned = (Q_centered @ R) + centroid_P
    return Q_aligned, R

def calculate_rmsd(P, Q_aligned):
    """محاسبه RMSD بین دو مجموعه نقطه همتراز شده"""
    diff = P - Q_aligned
    N = P.shape[0]
    rmsd = np.sqrt(np.sum(diff * diff) / N)
    return rmsd

def visualize_structures_complete(coords_A, coords_B, coords_B_aligned, rmsd, protein_names):
    """ویژوال کردن ساختارهای پروتئینی به صورت جداگانه و روی هم"""
    
    # ایجاد figure بزرگ با چند subplot
    fig = plt.figure(figsize=(20, 15))
    
    # تنظیمات پیش‌فرض subplot با مقادیر مورد نظر شما
    plt.subplots_adjust(
        left=0.045,
        bottom=0.055,
        right=0.955,
        top=0.936,
        wspace=0.212,
        hspace=0.221
    )
    
    # 1. نمای سه بعدی اصلی - روی هم
    ax1 = fig.add_subplot(331, projection='3d')
    ax1.plot(coords_A[:, 0], coords_A[:, 1], coords_A[:, 2], 
             'b-', linewidth=3, label=protein_names[0], alpha=0.8)
    ax1.scatter(coords_A[:, 0], coords_A[:, 1], coords_A[:, 2], 
                c='blue', s=20, alpha=0.6)
    
    ax1.plot(coords_B_aligned[:, 0], coords_B_aligned[:, 1], coords_B_aligned[:, 2], 
             'r-', linewidth=3, label=f'{protein_names[1]} (Aligned)', alpha=0.8)
    ax1.scatter(coords_B_aligned[:, 0], coords_B_aligned[:, 1], coords_B_aligned[:, 2], 
                c='red', s=20, alpha=0.6)
    
    ax1.set_xlabel('X (Å)')
    ax1.set_ylabel('Y (Å)')
    ax1.set_zlabel('Z (Å)')
    ax1.set_title(f'Overlay: {protein_names[0]} vs {protein_names[1]} | RMSD: {rmsd:.3f} Å')
    ax1.legend()
    ax1.view_init(elev=20, azim=45)
    
    
    # 2. نمای جانبی از overlay
    ax2 = fig.add_subplot(332, projection='3d')
    ax2.plot(coords_A[:, 0], coords_A[:, 1], coords_A[:, 2], 
             'b-', linewidth=2, label=protein_names[0], alpha=0.8)
    ax2.plot(coords_B_aligned[:, 0], coords_B_aligned[:, 1], coords_B_aligned[:, 2], 
             'r-', linewidth=2, label=f'{protein_names[1]} (Aligned)', alpha=0.8)
    
    # رسم خطوط بین اتم‌های متناظر
    for i in range(0, min(len(coords_A), len(coords_B_aligned)), 5):  # هر ۵ اتم یک خط
        ax2.plot([coords_A[i, 0], coords_B_aligned[i, 0]],
                [coords_A[i, 1], coords_B_aligned[i, 1]],
                [coords_A[i, 2], coords_B_aligned[i, 2]],
                'gray', alpha=0.3, linewidth=0.5)
    
    ax2.set_xlabel('X (Å)')
    ax2.set_ylabel('Y (Å)')
    ax2.set_zlabel('Z (Å)')
    ax2.set_title(f'Side View {protein_names[0]} vs {protein_names[1]}')
    ax2.legend()
    ax2.view_init(elev=10, azim=90)
    
    # 3. پروتئین اول به تنهایی
    ax3 = fig.add_subplot(333, projection='3d')
    ax3.plot(coords_A[:, 0], coords_A[:, 1], coords_A[:, 2], 
             'b-', linewidth=3, label=f'{protein_names[0]} (Original)', alpha=0.9)
    ax3.scatter(coords_A[:, 0], coords_A[:, 1], coords_A[:, 2], 
                c='blue', s=25, alpha=0.7)
    
    ax3.set_xlabel('X (Å)')
    ax3.set_ylabel('Y (Å)')
    ax3.set_zlabel('Z (Å)')
    ax3.set_title(f'{protein_names[0]} - Original Structure')
    ax3.legend()
    ax3.view_init(elev=20, azim=45)
    
    # 4. مقایسه ساختار اصلی پروتئین دوم با همتراز شده
    ax4 = fig.add_subplot(334, projection='3d')
    ax4.plot(coords_B[:, 0], coords_B[:, 1], coords_B[:, 2], 
             'green', linewidth=2, label=f'{protein_names[1]} (Original)', alpha=0.7)
    ax4.plot(coords_B_aligned[:, 0], coords_B_aligned[:, 1], coords_B_aligned[:, 2], 
             'red', linewidth=2, label=f'{protein_names[1]} (Aligned)', alpha=0.7)
    
    ax4.set_xlabel('X (Å)')
    ax4.set_ylabel('Y (Å)')
    ax4.set_zlabel('Z (Å)')
    ax4.set_title(f'{protein_names[1]} - Original vs Aligned')
    ax4.legend()
    ax4.view_init(elev=20, azim=45)
    
    # 5. پروتئین دوم به تنهایی (ساختار همتراز شده)
    ax5 = fig.add_subplot(335, projection='3d')
    ax5.plot(coords_B_aligned[:, 0], coords_B_aligned[:, 1], coords_B_aligned[:, 2], 
             'red', linewidth=3, label=f'{protein_names[1]} (Aligned)', alpha=0.9)
    ax5.scatter(coords_B_aligned[:, 0], coords_B_aligned[:, 1], coords_B_aligned[:, 2], 
                c='red', s=25, alpha=0.7)
    
    ax5.set_xlabel('X (Å)')
    ax5.set_ylabel('Y (Å)')
    ax5.set_zlabel('Z (Å)')
    ax5.set_title(f'{protein_names[1]} - Aligned Structure')
    ax5.legend()
    ax5.view_init(elev=20, azim=45)

    # 6. پروتئین دوم به تنهایی (ساختار اصلی)
    ax6 = fig.add_subplot(336, projection='3d')
    ax6.plot(coords_B[:, 0], coords_B[:, 1], coords_B[:, 2], 
             'green', linewidth=3, label=f'{protein_names[1]} (Original)', alpha=0.9)
    ax6.scatter(coords_B[:, 0], coords_B[:, 1], coords_B[:, 2], 
                c='green', s=25, alpha=0.7)
    
    ax6.set_xlabel('X (Å)')
    ax6.set_ylabel('Y (Å)')
    ax6.set_zlabel('Z (Å)')
    ax6.set_title(f'{protein_names[1]} - Original Structure')
    ax6.legend()
    ax6.view_init(elev=20, azim=45)

    # 7. هیستوگرام فاصله‌ها
    ax7 = fig.add_subplot(337)
    distances = np.linalg.norm(coords_A - coords_B_aligned, axis=1)
    n, bins, patches = ax7.hist(distances, bins=30, alpha=0.7, color='purple', edgecolor='black')
    ax7.set_xlabel('Distance (Å)')
    ax7.set_ylabel('Frequency')
    ax7.set_title('Distribution of Atomic Distances')
    ax7.grid(True, alpha=0.3)
    
    mean_distance = np.mean(distances)
    ax7.axvline(mean_distance, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_distance:.2f} Å')
    ax7.axvline(rmsd, color='blue', linestyle='--', linewidth=2, 
                label=f'RMSD: {rmsd:.2f} Å')
    ax7.legend()
    
    # 8. نمودار فاصله بر اساس شماره رزیدو
    ax8 = fig.add_subplot(338)
    residue_numbers = np.arange(1, len(distances) + 1)
    ax8.plot(residue_numbers, distances, 'orange', alpha=0.7, linewidth=1.5)
    ax8.fill_between(residue_numbers, distances, alpha=0.3, color='orange')
    ax8.set_xlabel('Residue Number')
    ax8.set_ylabel('Distance (Å)')
    ax8.set_title('Distance per Residue')
    ax8.grid(True, alpha=0.3)
    
    high_distance_threshold = np.percentile(distances, 75)
    high_dist_indices = distances > high_distance_threshold
    ax8.scatter(residue_numbers[high_dist_indices], distances[high_dist_indices], 
                color='red', s=30, zorder=5, label=f'Top 25% > {high_distance_threshold:.1f}Å')
    ax8.legend()
    
    # 9. نمودار تجمعی فاصله
    ax9 = fig.add_subplot(339)
    sorted_distances = np.sort(distances)
    cumulative_prob = np.arange(1, len(sorted_distances) + 1) / len(sorted_distances)
    ax9.plot(sorted_distances, cumulative_prob, 'brown', linewidth=2)
    ax9.set_xlabel('Distance (Å)')
    ax9.set_ylabel('Cumulative Probability')
    ax9.set_title('Cumulative Distance Distribution')
    ax9.grid(True, alpha=0.3)
    
    ax9.axhline(0.5, color='red', linestyle='--', alpha=0.7, label='Median')
    ax9.axvline(np.median(distances), color='red', linestyle='--', alpha=0.7)
    ax9.legend()
    
    plt.show()
    
    return distances

# برنامه اصلی
if __name__ == "__main__":
    # مسیرهای مستقیم فایل‌های PDB
    # 5TIM - 1TIM
    # 1A01 - 1A00
    pdb_file_1 = r"C:\Users\ARNZ\Desktop\1A00.pdb"
    pdb_file_2 = r"C:\Users\ARNZ\Desktop\1A01.pdb"
    
    protein_names = ["1A00", "1A01"]
    
    print("=== Protein 3D Structure Comparison ===")
    print(f"Protein 1: {protein_names[0]}")
    print(f"Protein 2: {protein_names[1]}")
    print("-" * 50)
    
    # بررسی وجود فایل‌ها
    if not os.path.exists(pdb_file_1):
        print(f"ERROR: File {pdb_file_1} does not exist!")
        exit()
    if not os.path.exists(pdb_file_2):
        print(f"ERROR: File {pdb_file_2} does not exist!")
        exit()
    
    # خواندن ساختارها
    print("Reading PDB files...")
    coords_protein_A = read_pdb(pdb_file_1)
    coords_protein_B_original = read_pdb(pdb_file_2)  # ذخیره ساختار اصلی
    
    if coords_protein_A is None or coords_protein_B_original is None:
        print("Failed to read PDB files. Exiting...")
        exit()
    
    print(f"Protein A: {len(coords_protein_A)} CA atoms found")
    print(f"Protein B: {len(coords_protein_B_original)} CA atoms found")
    
    # کپی از مختصات اصلی برای نمایش جداگانه
    coords_protein_B = coords_protein_B_original.copy()
    
    # بررسی طول پروتئین‌ها
    if len(coords_protein_A) != len(coords_protein_B):
        print(f"WARNING: Proteins have different lengths ({len(coords_protein_A)} vs {len(coords_protein_B)})")
        print("Using only the first N atoms where N is the minimum length")
        min_len = min(len(coords_protein_A), len(coords_protein_B))
        coords_protein_A = coords_protein_A[:min_len]
        coords_protein_B = coords_protein_B[:min_len]
        print(f"Using {min_len} atoms for comparison")
    
    # همترازی ساختاری
    print("Performing structural alignment...")
    coords_protein_B_aligned, rotation_matrix = kabsch_algorithm(coords_protein_A, coords_protein_B)
    
    # محاسبه RMSD
    print("Calculating RMSD...")
    rmsd_value = calculate_rmsd(coords_protein_A, coords_protein_B_aligned)
    
    print("-" * 50)
    print("RESULTS:")
    print(f"RMSD between structures: {rmsd_value:.3f} Angstroms")
    print(f"Number of residues compared: {len(coords_protein_A)}")
    
    # تفسیر نتیجه
    print("\nRMSD Interpretation:")
    if rmsd_value < 1.0:
        print("✓ Very high structural similarity - Likely same fold")
    elif rmsd_value < 2.0:
        print("✓ High structural similarity - Very similar structures")
    elif rmsd_value < 3.0:
        print("○ Moderate structural similarity - Related structures")
    elif rmsd_value < 4.0:
        print("△ Low structural similarity - Some common features")
    else:
        print("✗ Very low structural similarity - Different folds")
    
    print("-" * 50)
    print("Generating complete visualization...")
    
    # ویژوال‌سازی کامل با نمایش جداگانه
    distances = visualize_structures_complete(
        coords_protein_A, 
        coords_protein_B_original[:len(coords_protein_A)] if len(coords_protein_B_original) >= len(coords_protein_A) else coords_protein_B_original,
        coords_protein_B_aligned, 
        rmsd_value, 
        protein_names
    )
    
    # اطلاعات آماری اضافی
    print("\nAdditional Statistics:")
    print(f"Mean distance: {np.mean(distances):.3f} Å")
    print(f"Median distance: {np.median(distances):.3f} Å")
    print(f"Max distance: {np.max(distances):.3f} Å")
    print(f"Min distance: {np.min(distances):.3f} Å")
    print(f"Standard deviation: {np.std(distances):.3f} Å")
    print(f"25th percentile: {np.percentile(distances, 25):.3f} Å")
    print(f"75th percentile: {np.percentile(distances, 75):.3f} Å")