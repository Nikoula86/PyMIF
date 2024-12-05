import glob
import subprocess

scripts = glob.glob("test_*.py")
scripts = [s for s in scripts if not "_all.py" in s]
print(scripts)

for i, script in enumerate(scripts):
    print(f"### {i}. Running {script}")
    subprocess.run(["python", script])
