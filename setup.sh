#!/usr/bin/env bash
set -e

echo "======================================"
echo " Quantum Bio Recovery Engine — Setup "
echo "======================================"

# Check Python
if ! command -v python3 &> /dev/null; then
  echo "❌ Python3 not found. Please install Python 3.10+"
  exit 1
fi

PY_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
echo "🐍 Python version: $PY_VERSION"

# Create virtual environment
if [ ! -d "qbio-env" ]; then
  echo "📦 Creating virtual environment: qbio-env"
  python3 -m venv qbio-env
else
  echo "📦 Virtual environment already exists: qbio-env"
fi

# Activate environment
echo "⚡ Activating virtual environment"
source qbio-env/bin/activate

# Upgrade pip
echo "⬆️  Upgrading pip"
pip install --upgrade pip

# Install dependencies
echo "📥 Installing dependencies from requirements.txt"
pip install -r requirements.txt

echo ""
echo "✅ Setup complete!"
echo ""
echo "Next steps:"
echo "  source qbio-env/bin/activate"
echo "  python -m experiments.benchmarks.scale_benchmark"
echo "  python experiments/benchmarks/plot_scaling.py"
echo ""
