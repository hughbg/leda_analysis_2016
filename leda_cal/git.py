import os
import subprocess

def get_hash():
    """
    Find and return the hash of the last git commit.  Raises a RuntimeError if
    there is a problem determing the hash.
    """
    
    p = subprocess.Popen(['git', 'rev-parse', '--verify', 'HEAD'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.wait() != 0:
        raise RuntimeError(error)
        
    return output.strip()


def get_branch():
    """
    Find and return the name of the current git branch.  Raises a RuntimeError
    if there is a problem determing the hash.
    """
    
    p = subprocess.Popen(['git', 'rev-parse', '--abbrev-ref', 'HEAD'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.wait() != 0:
        raise RuntimeError(error)
        
    return output.strip()


def is_branch_dirty(exclude_exts=None):
    """
    Find out if there are any unstagged commit in the current git branch.  
    Raises a RuntimeError if there is a problem determing the state of the
    branch.
    
    The 'exclude_exts' keyword is used to filter out unstagged commits on 
    files based on their extensions
    """
    
    p = subprocess.Popen(['git', 'diff-index', 'HEAD', '--'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = p.communicate()
    if p.wait() != 0:
        raise RuntimeError(error)
        
    # Check the output for details
    if exclude_exts is not None:
        if type(exclude_exts) is not list:
            exclude_exts = [exclude_exts,]
            
        dirty = False
        output = output.split('\n')[:-1]
        for line in output:
            filename = line.split()[-1]
            _, ext = os.path.splitext(filename)
            if ext not in exclude_exts:
                dirty = True
    else:
        dirty = True if len(output) > 3 else False
        
    return dirty


def get_repo_fingerprint():
    """
    Return a string that encodes the current state of the git repository so 
    that files created by the analysis can be tagged with the software used.
    Raises a RuntimeError if there is a problem determing the state of the repo.
    """
    
    branch = get_branch()
    ghash = get_hash()
    dirty = '_dirty' if is_branch_dirty(exclude_exts=['.pdf', '.png']) else ''
    
    return "%s%s_%s" % (branch, dirty, ghash)


if __name__ == "__main__":
    print get_branch()
    print get_hash()
    print is_branch_dirty()
    print get_repo_fingerprint()