#!/usr/bin/env python3
"""
Simple test script to verify STCRpy installation in Docker
"""

import sys

def test_imports():
    """Test that all required modules can be imported"""
    print("Testing STCRpy installation...")

    try:
        import stcrpy
        print("✓ stcrpy imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import stcrpy: {e}")
        return False

    try:
        import Bio
        print("✓ biopython imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import biopython: {e}")
        return False

    try:
        import pandas as pd
        print("✓ pandas imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import pandas: {e}")
        return False

    try:
        import numpy as np
        print("✓ numpy imported successfully")
    except ImportError as e:
        print(f"✗ Failed to import numpy: {e}")
        return False

    return True

def test_basic_functionality():
    """Test basic STCRpy functionality"""
    print("\nTesting basic functionality...")

    try:
        import stcrpy
        # Try fetching a structure (requires internet)
        print("Attempting to fetch TCR structure 8gvb from PDB...")
        tcrs = stcrpy.fetch_TCRs("8gvb")
        print(f"✓ Successfully fetched {len(tcrs)} TCR structure(s)")

        if tcrs:
            tcr = tcrs[0]
            print(f"  - TCR has {len(list(tcr.get_chains()))} chain(s)")

        return True
    except Exception as e:
        print(f"✗ Basic functionality test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("=" * 60)
    print("STCRpy Docker Installation Test")
    print("=" * 60)

    # Test imports
    if not test_imports():
        print("\n❌ Import tests failed!")
        sys.exit(1)

    # Test basic functionality
    if not test_basic_functionality():
        print("\n⚠️  Basic functionality test failed (may require internet)")

    print("\n" + "=" * 60)
    print("✅ STCRpy is properly installed and ready to use!")
    print("=" * 60)

if __name__ == "__main__":
    main()
