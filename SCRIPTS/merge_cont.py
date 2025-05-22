import os
import sys
from PIL import Image

def merge_images(png1_path, png2_path, output_path):
    image1 = Image.open(png1_path)
    image2 = Image.open(png2_path)
    width1, height1 = image1.size
    width2, height2 = image2.size
    total_width = width1 + width2
    total_height = max(height1, height2)
    merged_image = Image.new('RGB', (total_width, total_height), color='white')
    merged_image.paste(image1, (0, 0))
    merged_image.paste(image2, (width1, 0))
    merged_image.save(output_path)

def merge_plots(directory):
    file_dict = {}
    for file in os.listdir(directory):
        if file.endswith(".png"):
             parts = file.split('_')
             if len(parts) >= 3:
                key = '_'.join(parts[:2])
                if key not in file_dict:
                    file_dict[key] = [file]
                else:
                    file_dict[key].append(file)

    # Traverse the dictionary and merge the corresponding files
    for key, value in file_dict.items():
        if len(value) == 2:
            # Separate files by type
            plot_phase_file = next(f for f in value if 'plot_phase' in f)
            plot_file = next(f for f in value if 'plot' in f and 'plot_phase' not in f)
            output_path = os.path.join(directory, f"{key}_merged.png")
            merge_images(os.path.join(directory, plot_file), os.path.join(directory, plot_phase_file), output_path)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py directory_path")
        sys.exit(1)
        
    directory = sys.argv[1]
    merge_plots(directory)
