import numpy as np
import os
import os.path as op
import sys

rl = sys.argv[1]
PRD = os.environ['PRD']
FS = os.environ['FS']
SUBJ_ID = os.environ['SUBJ_ID']

def read_annot(fname):

    """Read a Freesurfer annotation from a .annot file.
    Note : Copied from nibabel
    Parameters
    ----------
    fname : str
        Path to annotation file
    Returns
    -------
    annot : numpy array, shape=(n_verts)
        Annotation id at each vertex
    ctab : numpy array, shape=(n_entries, 5)
        RGBA + label id colortable array
    names : list of str
        List of region names as stored in the annot file
    """
    if not op.isfile(fname):
        dir_name = op.split(fname)[0]
        if not op.isdir(dir_name):
            raise IOError('Directory for annotation does not exist: %s',
                          fname)
        cands = os.listdir(dir_name)
        cands = [c for c in cands if '.annot' in c]
        if len(cands) == 0:
            raise IOError('No such file %s, no candidate parcellations '
                          'found in directory' % fname)
        else:
            raise IOError('No such file %s, candidate parcellations in '
                          'that directory: %s' % (fname, ', '.join(cands)))
    with open(fname, "rb") as fid:
        n_verts = np.fromfile(fid, '>i4', 1)[0]
        data = np.fromfile(fid, '>i4', n_verts * 2).reshape(n_verts, 2)
        annot = data[data[:, 0], 1]
        ctab_exists = np.fromfile(fid, '>i4', 1)[0]
        if not ctab_exists:
            raise Exception('Color table not found in annotation file')
        n_entries = np.fromfile(fid, '>i4', 1)[0]
        if n_entries > 0:
            length = np.fromfile(fid, '>i4', 1)[0]
            orig_tab = np.fromfile(fid, '>c', length)
            orig_tab = orig_tab[:-1]
            names = list()
            ctab = np.zeros((n_entries, 5), np.int)
            for i in range(n_entries):
                name_length = np.fromfile(fid, '>i4', 1)[0]
                name = np.fromfile(fid, "|S%d" % name_length, 1)[0]
                names.append(name)
                ctab[i, :4] = np.fromfile(fid, '>i4', 4)
                ctab[i, 4] = (ctab[i, 0] + ctab[i, 1] * (2 ** 8) +
                              ctab[i, 2] * (2 ** 16) +
                              ctab[i, 3] * (2 ** 24))
        else:
            ctab_version = -n_entries
            if ctab_version != 2:
                raise Exception('Color table version not supported')
            n_entries = np.fromfile(fid, '>i4', 1)[0]
            ctab = np.zeros((n_entries, 5), np.int)
            length = np.fromfile(fid, '>i4', 1)[0]
            np.fromfile(fid, "|S%d" % length, 1)  # Orig table path
            entries_to_read = np.fromfile(fid, '>i4', 1)[0]
            names = list()
            for i in range(entries_to_read):
                np.fromfile(fid, '>i4', 1)  # Structure
                name_length = np.fromfile(fid, '>i4', 1)[0]
                name = np.fromfile(fid, "|S%d" % name_length, 1)[0]
                names.append(name)
                ctab[i, :4] = np.fromfile(fid, '>i4', 4)
                ctab[i, 4] = (ctab[i, 0] + ctab[i, 1] * (2 ** 8) +
                              ctab[i, 2] * (2 ** 16))
        # convert to more common alpha value
        ctab[:, 3] = 255 - ctab[:, 3]
    return annot, ctab, names


L, _, _ = read_annot(os.path.join(FS, SUBJ_ID, 'label', rl + '.aparc.annot'))
ref_table = np.loadtxt(os.path.join('share', rl + '_ref_table.txt'))
os.chdir(os.path.join(PRD, 'surface'))
vl = np.loadtxt(rl + '_vertices_low.txt')  # vertices low
vh = np.loadtxt(rl + '_vertices_high.txt')  # vertices high
reg_map = []
for vli in vl:
    pos = np.argmin(np.sum(np.abs(vh - vli), 1))
    find_tab = np.nonzero(ref_table[:, 5] == L[pos])[0][0]
    reg_map.append(ref_table[find_tab, 4])
np.savetxt(rl + '_region_mapping_low_not_corrected.txt', reg_map)
