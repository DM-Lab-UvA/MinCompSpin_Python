import sys
import os

# Add the packages directory to the path
package_dir = os.path.join(os.path.dirname(__file__), 'packages')
if package_dir not in sys.path:
    print(package_dir)
    sys.path.append(package_dir)

# Now import your module
try:
    import MinCompSpin
    print("Module imported successfully!")
except ImportError as e:
    print(f"Failed to import module: {e}")
