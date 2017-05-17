"""
Adds main qcore directory to PATH.
Add QuakeCoRE Library to PATH
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
