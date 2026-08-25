import json
import os

BASE_DIR = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)

TARGET_DIR = os.path.join(BASE_DIR, "data", "targets")

print("🧪 TARGET_DIR =", TARGET_DIR)
print("🧪 Files in TARGET_DIR =", os.listdir(TARGET_DIR) if os.path.exists(TARGET_DIR) else "DIR NOT FOUND")


def load_targets(disease):
    path = os.path.join(TARGET_DIR, f"{disease}.json")
    print("🧪 Looking for:", path)

    try:
        with open(path, "r") as f:
            return json.load(f).get("targets", [])
    except Exception as e:
        print("🧪 Load error:", e)
        return []
