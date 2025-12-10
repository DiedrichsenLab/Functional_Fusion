import os
import subprocess

new_directory = 'Y:/data/ExternalOpenData/HCP_UR100_new/tasktest'
if not os.path.exists(new_directory):
    new_directory = '/cifs/diedrichsen/data/ExternalOpenData/HCP_UR100_new/tasktest'

def run_all_glms(directory):
    """
    Finds every .fsf file inside the directory, updates its outputdir,
    and runs FEAT on it.
    """
    for root, dirs, files in os.walk(directory):
        for file in files:
            if not file.endswith(".fsf"):
                continue

            fsf_path = os.path.join(root, file)

            # Build updated output directory
            fsf_dir = os.path.dirname(fsf_path)
            base = os.path.splitext(file)[0]
            output_dir = os.path.join(fsf_dir, f"{base}.feat")

            # ---- Update outputdir in .fsf ----
            updated_lines = []
            with open(fsf_path, "r") as f:
                for line in f:
                    if line.startswith("set fmri(outputdir)"):
                        updated_lines.append(f'set fmri(outputdir) "{output_dir}"\n')
                    else:
                        updated_lines.append(line)

            with open(fsf_path, "w") as f:
                f.writelines(updated_lines)

            print(f"[OK] Updated outputdir → {output_dir}")

            # ---- Run FEAT ----
            try:
                print(f"[RUN] feat {fsf_path}")
                subprocess.run(["feat", fsf_path], check=True)
                print(f"[DONE] FEAT complete → {fsf_path}")

            except subprocess.CalledProcessError as e:
                print(f"[ERROR] FEAT failed on {fsf_path}: {e}")
            except Exception as e:
                print(f"[ERROR] Unexpected error: {e}")

    print("\n=== All GLMs finished ===")



if __name__ == "__main__":
    run_all_glms(new_directory)



    

    