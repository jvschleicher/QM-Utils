#!/usr/bin/env python3
"""
Script to efficiently restart failed ORCA jobs
"""

import os
import glob
import re
from pathlib import Path
import argparse


def check_job_status(outfile):
    """
    Check if ORCA job completed successfully or failed.
    
    Returns:
        tuple: (status, failure_type)
            status: 'success', 'failed', 'running', or 'not_started'
            failure_type: None, 'scf', 'opt', 'timeout', 'other'
    """
    if not os.path.exists(outfile):
        return 'not_started', None
    
    try:
        with open(outfile, 'r') as f:
            # Read last 10 lines to check for termination
            lines = f.readlines()
            if len(lines) < 5:
                return 'not_started', None
            
            last_lines = ''.join(lines[-10:])
            full_content = ''.join(lines)
    except:
        return 'not_started', None
    
    # Check for successful completion (must be in last few lines)
    if 'ORCA TERMINATED NORMALLY' in last_lines:
        return 'success', None
    
    # If file has content but no normal termination, it failed
    if len(full_content) < 100:
        return 'not_started', None
    
    # Check for various failure modes to categorize
    failure_type = 'other'
    if 'SCF NOT CONVERGED' in full_content or 'SCF CONVERGENCE ... NOT CONVERGED' in full_content:
        failure_type = 'scf'
    elif 'THE OPTIMIZATION DID NOT CONVERGE' in full_content:
        failure_type = 'opt'
    elif 'ABORTING THE RUN' in full_content:
        failure_type = 'other'
    
    return 'failed', failure_type


def read_input_file(inpfile):
    """Read and parse ORCA input file."""
    with open(inpfile, 'r') as f:
        return f.read()


def extract_input_line_from_output(outfile):
    """
    Extract the original input line from ORCA output file.
    Looks for the INPUT FILE section.
    
    Returns:
        str: Input line starting with '!' or None
    """
    if not os.path.exists(outfile):
        return None
    
    with open(outfile, 'r') as f:
        lines = f.readlines()
    
    in_input_section = False
    lines_read_in_section = 0
    for i, line in enumerate(lines):
        if 'INPUT FILE' in line:
            in_input_section = True
            continue
        
        if in_input_section:
            lines_read_in_section += 1
            
            # Look for lines like "|  1> ! HF-3c opt numfreq TightSCF"
            if '|' in line and '>' in line and '!' in line:
                # Extract the part after '>'
                parts = line.split('>', 1)
                if len(parts) == 2:
                    input_line = parts[1].strip()
                    if input_line.startswith('!'):
                        return input_line
            
            # Only stop at section boundary after we've read at least a few lines
            # (to skip the ======== line right after INPUT FILE)
            if lines_read_in_section > 3:
                if line.strip() and line.strip()[0] == '=' and line.strip()[-1] == '=':
                    break
    
    return None


def generate_restart_input(inpfile, outfile, failure_type, custom_input_line=None):
    """
    Generate restart input file with modified settings based on failure type.
    
    Args:
        inpfile: Original input file path
        outfile: Output file path to extract geometry
        failure_type: Type of failure ('scf', 'opt', etc.)
        custom_input_line: Custom input line to use (if None, extract from .out file)
        
    Returns:
        str: Content for restart input file or None on error
    """
    try:
        # Get input line - use custom if provided, otherwise extract from .out file
        if custom_input_line:
            settings = custom_input_line.strip()
        else:
            settings = extract_input_line_from_output(outfile)
            if not settings:
                # Fallback to reading from .inp file if extraction fails
                if os.path.exists(inpfile):
                    original_content = read_input_file(inpfile)
                    lines = original_content.split('\n')
                    for line in lines:
                        if line.strip().startswith('!'):
                            settings = line.strip()
                            break
                if not settings:
                    print(f"\n    ERROR: Could not extract input line from {outfile}")
                    return None
        
        # Extract charge and multiplicity from output or input file
        charge = 0
        mult = 1
        
        # First try to extract from output file
        with open(outfile, 'r') as f:
            for line in f:
                # Look for the xyz line in INPUT FILE section
                if '* xyz' in line.lower() and '|' in line:
                    # Format: "|  2> * xyz 0 1"
                    parts = line.split('>')
                    if len(parts) > 1:
                        xyz_parts = parts[1].split()
                        if len(xyz_parts) >= 4:
                            charge = xyz_parts[2]
                            mult = xyz_parts[3]
                            break
        
        # Fallback: try to get from input file if it exists and we didn't find it
        if (charge == 0 and mult == 1) and os.path.exists(inpfile):
            original_content = read_input_file(inpfile)
            lines = original_content.split('\n')
            for line in lines:
                if '* xyz' in line.lower():
                    xyz_line_parts = line.split()
                    if len(xyz_line_parts) >= 4:
                        charge = xyz_line_parts[2]
                        mult = xyz_line_parts[3]
                        break
        
        # Get last geometry from output file if available
        last_geometry = extract_last_geometry(outfile)
        
        if not last_geometry:
            print(f"\n    ERROR: Could not extract geometry from {outfile}")
            return None
        
        if failure_type == 'scf':
            # Add more robust SCF convergence settings
            if 'SlowConv' not in settings:
                settings = settings.replace('!', '! SlowConv', 1)
            if 'VeryTightSCF' in settings:
                settings = settings.replace('VeryTightSCF', 'TightSCF')
            if 'TightSCF' not in settings and 'SlowConv' not in settings:
                settings = settings.replace('!', '! SlowConv', 1)
                
        elif failure_type == 'opt':
            # Relax optimization convergence
            if 'TightOpt' in settings:
                settings = settings.replace('TightOpt', 'NormalOpt')
            if 'VeryTightOpt' in settings:
                settings = settings.replace('VeryTightOpt', 'NormalOpt')
        
        # Build restart input
        restart_content = f"{settings}\n\n"
        restart_content += f"* xyz {charge} {mult}\n"
        restart_content += last_geometry
        restart_content += "*\n"
        
        return restart_content
        
    except Exception as e:
        print(f"\n    ERROR generating restart for {outfile}: {str(e)}")
        return None


def extract_last_geometry(outfile):
    """
    Extract the last Cartesian geometry from ORCA output file.
    Uses the same parsing logic as parse_orca_xyz.py
    
    Returns:
        str: Formatted geometry or None
    """
    if not os.path.exists(outfile):
        return None
    
    frames = []
    
    with open(outfile, 'r') as f:
        lines = iter(f)
        for line in lines:
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                current = []
                next(lines)  # Skip dashed line
                # Process atoms immediately after the dashed line
                for atom_line in lines:
                    if atom_line.strip() == "":
                        if current:
                            frames.append(current)
                        current = None
                        break
                    
                    parts = atom_line.split()
                    if len(parts) == 4:
                        current.append(parts)
    
    if not frames:
        return None
    
    # Format the last geometry
    last_frame = frames[-1]
    geometry_lines = []
    for el, x, y, z in last_frame:
        geometry_lines.append(f"{el:2s} {x:>12s} {y:>12s} {z:>12s}")
    
    return '\n'.join(geometry_lines) + '\n'


def process_directory(directory, pattern="*.out", dry_run=True, suffix="_restart", custom_input_line=None):
    """
    Process all ORCA output files in a directory and generate restart inputs.
    
    Args:
        directory: Directory to scan
        pattern: Pattern to match output files
        dry_run: If True, only report status without creating restart files
        suffix: Suffix for restart input files
        custom_input_line: Custom input line to use for all restarts (if None, extract from each .out file)
    """
    out_files = glob.glob(os.path.join(directory, pattern))
    
    print(f"Scanning {directory} for ORCA jobs...")
    print(f"Found {len(out_files)} output files\n")
    
    summary = {'success': 0, 'failed': 0, 'running': 0, 'not_started': 0}
    failures = {'scf': 0, 'opt': 0, 'other': 0}
    restart_files = []
    
    for outfile in out_files:
        status, failure_type = check_job_status(outfile)
        summary[status] += 1
        
        basename = os.path.splitext(os.path.basename(outfile))[0]
        inpfile = os.path.join(directory, f"{basename}.inp")
        
        print(f"{basename:50s} {status:12s}", end='')
        
        if status == 'failed':
            print(f" ({failure_type})", end='')
            failures[failure_type] += 1
            
            if not dry_run:
                restart_content = generate_restart_input(inpfile, outfile, failure_type, custom_input_line)
                if restart_content:
                    restart_file = os.path.join(directory, f"{basename}{suffix}.inp")
                    with open(restart_file, 'w') as f:
                        f.write(restart_content)
                    restart_files.append(restart_file)
                    print(f" -> {basename}{suffix}.inp", end='')
        
        print()
    
    print("\n" + "="*80)
    print("SUMMARY:")
    print(f"  Success: {summary['success']}")
    print(f"  Failed:  {summary['failed']}")
    if failures['scf'] > 0:
        print(f"    - SCF convergence: {failures['scf']}")
    if failures['opt'] > 0:
        print(f"    - Optimization: {failures['opt']}")
    if failures['other'] > 0:
        print(f"    - Other: {failures['other']}")
    print(f"  Running: {summary['running']}")
    print(f"  Not started: {summary['not_started']}")
    
    if not dry_run and restart_files:
        print(f"\nCreated {len(restart_files)} restart input files")
    elif dry_run and summary['failed'] > 0:
        print("\nRun with --execute to create restart input files")
    
    return restart_files


def main():
    parser = argparse.ArgumentParser(
        description='Efficiently restart failed ORCA jobs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan directory and report status (dry run)
  python restart_orca_jobs.py /path/to/jobs
  
  # Create restart input files for failed jobs
  python restart_orca_jobs.py /path/to/jobs --execute
  
  # Custom output pattern and restart suffix
  python restart_orca_jobs.py /path/to/jobs --pattern "calc_*.out" --suffix "_v2"
        """
    )
    
    parser.add_argument('directory', help='Directory containing ORCA output files')
    parser.add_argument('--pattern', default='*.out', help='Pattern to match output files (default: *.out)')
    parser.add_argument('--suffix', default='_restart', help='Suffix for restart input files (default: _restart)')
    parser.add_argument('--execute', action='store_true', help='Create restart input files (default: dry run)')
    parser.add_argument('--input-line', dest='input_line', help='Custom input line for all restarts (e.g., "! HF-3c opt"). If not provided, extracts from each .out file')
    
    args = parser.parse_args()
    
    if not os.path.isdir(args.directory):
        print(f"Error: {args.directory} is not a valid directory")
        return 1
    
    process_directory(
        args.directory,
        pattern=args.pattern,
        dry_run=not args.execute,
        suffix=args.suffix,
        custom_input_line=args.input_line
    )
    
    return 0


if __name__ == "__main__":
    exit(main())
