#!/usr/bin/env python

strand_resids = []
strand_beads = []

def gen_pdb_file():
    '''generate CafeMol readable pdb files from atomistic.pdb'''
    all_resids = []
    strand_heads = []
    for strnd in strand_resids:
        all_resids.extend(strnd)
        strand_heads.append(strnd[0])

    out_pdb_file = open("../dna.pdb", 'w')
    with open("atomistic.pdb", 'r') as fin:
        for lines in fin:
            if lines.startswith('ATOM  '):
                resid = int(lines[22:26])
                atom_name = lines[12:16].strip()
                if resid not in all_resids:
                    continue
                if resid in strand_heads and atom_name in ['P', 'OP1', 'OP2']:
                    if resid != strand_heads[0] and atom_name == 'P':
                        out_pdb_file.write('TER\n')
                    continue
                out_pdb_file.write(lines)
        out_pdb_file.write('END\n')
    out_pdb_file.close()


def gen_ninfo_file():
    # ==================== generate resid_num trans dict ====================
    all_beads = []
    bead_to_strand = {}
    for i, strnd in enumerate(strand_beads):
        all_beads.extend(strnd)
        for j in strnd:
            bead_to_strand[j] = i
    # beads_renum = {j:i+1 for i, j in enumerate(all_beads)}
    beads_renum = {}
    for i, j in enumerate(all_beads):
        beads_renum[j] = i + 1

    # ==================== read parameters from list files ====================
    bond_dict = {}
    angl_dict = {}
    dihe_dict = {}

    f_in_bond = open("in00_bond.list", 'r')
    f_in_angl = open("in00_angl.list", 'r')
    f_in_dihe = open("in00_dihe.list", 'r')

    for lines in f_in_bond:
        words = lines.split()
        imp1, imp2 = int(words[0]), int(words[1])
        if imp1 not in all_beads or imp2 not in all_beads:
            continue
        bd_nat = float(words[2])
        coef_bd = float(words[3])
        bond_dict[(imp1, imp2)] = [bd_nat, coef_bd]

    for lines in f_in_angl:
        words = lines.split()
        imp1, imp2, imp3 = int(words[0]), int(words[1]), int(words[2])
        if imp1 not in all_beads or imp2 not in all_beads or imp3 not in all_beads:
            continue
        ba_nat = float(words[3])
        coef_ba = float(words[4])
        angl_dict[(imp1, imp2, imp3)] = [ba_nat, coef_ba]

    for lines in f_in_dihe:
        words = lines.split()
        imp1, imp2, imp3, imp4 = int(words[0]), int(words[1]), int(words[2]), int(words[3])
        if imp1 not in all_beads or imp2 not in all_beads:
            continue
        if imp3 not in all_beads or imp4 not in all_beads:
            continue
        dih_nat = float(words[4])
        coef_prd = float(words[5])
        coef_gau = float(words[6])
        coef_sig = float(words[7])
        coef_type = int(words[8])
        dihe_dict[(imp1, imp2, imp3, imp4)] = [dih_nat, coef_prd, coef_gau, coef_sig, coef_type]

    f_in_bond.close()
    f_in_angl.close()
    f_in_dihe.close()


    # ==================== output params to ninfo files ====================
    # ---------- output strings ----------
    bond_string = "bond   {0:4d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d} {6:6d} {7:12.4f} {8:12.4f} {9:12.4f} {10:12.4f}\n"
    angl_string = "angl {0:6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d} {6:6d} {7:6d} {8:6d} {9:12.4f} {10:12.4f} {11:12.4f} {12:12.4f}\n"
    dihe_string = "dihd {0:6d} {1:6d} {2:6d} {3:6d} {4:6d} {5:6d} {6:6d} {7:6d} {8:6d} {9:6d} {10:6d} {11:12.4f} {12:12.4f} {13:12.4f} {14:12.4f} {15:12.4f} N2P{16:1d}\n"

    # ---------- out files ----------
    f_info_file = []
    for i in range(len(strand_resids)):
        f_info_file_name = "../strand" + str(i + 1) + ".ninfo"
        f_info_file.append(open(f_info_file_name, 'w'))

    # -------------------- BOND --------------------
    for f in f_info_file:
        f.write("<<<< native bond length\n")

    i_count = 0
    bond_keys = sorted(bond_dict.keys())
    for k in bond_keys:
        imp1, imp2 = k[0], k[1]
        snd1, snd2 = bead_to_strand[imp1], bead_to_strand[imp2]
        if snd1 != snd2:
            print("ERROR in bond params:", imp1, ' ', imp2)
            continue
        i_count += 1
        v = bond_dict[k]
        new_imp1, new_imp2 = beads_renum[imp1], beads_renum[imp2]
        loc_imp1 = strand_beads[snd1].index(imp1) + 1
        loc_imp2 = strand_beads[snd2].index(imp2) + 1
        f_info_file[snd1].write(bond_string.format(i_count, snd1+1, snd2+1, \
                                                   new_imp1, new_imp2, \
                                                   loc_imp1, loc_imp2, \
                                                   v[0], 1.0, 1.0, v[1]))
        # f_info_file[snd1].write("test")
    for f in f_info_file:
        f.write(">>>>\n\n")

    # -------------------- ANGLE --------------------
    for f in f_info_file:
        f.write("<<<< native bond angles\n")

    i_count = 0
    angl_keys = sorted(angl_dict.keys())
    for k in angl_keys:
        imp1, imp2, imp3 = k[0], k[1], k[2]
        snd1, snd2, snd3 = bead_to_strand[imp1], bead_to_strand[imp2], bead_to_strand[imp3]
        if not snd1 == snd2 == snd3:
            print("ERROR in angle params:", imp1, ' ', imp2, ' ', imp3)
            continue
        i_count += 1
        v = angl_dict[k]
        new_imp1, new_imp2, new_imp3 = beads_renum[imp1], beads_renum[imp2], beads_renum[imp3]
        loc_imp1 = strand_beads[snd1].index(imp1) + 1
        loc_imp2 = strand_beads[snd2].index(imp2) + 1
        loc_imp3 = strand_beads[snd3].index(imp3) + 1
        f_info_file[snd1].write(angl_string.format(i_count, snd1+1, snd1+1, \
                                                   new_imp1, new_imp2, new_imp3, \
                                                   loc_imp1, loc_imp2, loc_imp3, \
                                                   v[0], 1.0, 1.0, v[1]))
        # f_info_file[snd1].write("test")
    for f in f_info_file:
        f.write(">>>>\n\n")

    # -------------------- DIHEDRAL ANGLE --------------------
    for f in f_info_file:
        f.write("<<<< native dihedral angles\n")

    i_count = 0
    dihe_keys = sorted(dihe_dict.keys())
    for k in dihe_keys:
        imp1, imp2, imp3, imp4 = k[0], k[1], k[2], k[3]
        snd1, snd2 = bead_to_strand[imp1], bead_to_strand[imp2]
        snd3, snd4 = bead_to_strand[imp3], bead_to_strand[imp4]
        if not snd1 == snd2 == snd3 == snd4:
            print("ERROR in dihedral angle params:", imp1, ' ', imp2, ' ', imp3, ' ', imp4)
            continue
        i_count += 1
        v = dihe_dict[k]
        new_imp1, new_imp2 = beads_renum[imp1], beads_renum[imp2]
        new_imp3, new_imp4 = beads_renum[imp3], beads_renum[imp4]
        loc_imp1 = strand_beads[snd1].index(imp1) + 1
        loc_imp2 = strand_beads[snd2].index(imp2) + 1
        loc_imp3 = strand_beads[snd3].index(imp3) + 1
        loc_imp4 = strand_beads[snd4].index(imp4) + 1
        f_info_file[snd1].write(dihe_string.format(i_count, snd1+1, snd1+1, \
                                                   new_imp1, new_imp2, \
                                                   new_imp3, new_imp4, \
                                                   loc_imp1, loc_imp2, \
                                                   loc_imp3, loc_imp4, \
                                                   v[0], 1.0, 1.0, v[2], v[3], v[4]))
    for f in f_info_file:
        f.write(">>>>\n\n")
        f.close()

if __name__ == '__main__':
    # read in params from nickgap.inp;
    # nickgap.inp format:
    # ----------------------------------------
    # BDNALEN 10
    # STRAND1 1 to 5
    # STRAND2 6 to 10
    # STRAND3 11 to 20
    # ----------------------------------------

    import sys
    strand_info_fname = sys.argv[1]
    with open(strand_info_fname, 'r') as f:
        n_bp = 0                # number of nucleotides in one strand
        for lines in f:
            if lines.startswith('BDNALEN'):
                n_bp = int(lines.split()[1])
            if lines.startswith('STRAND'):
                words = lines.split()
                st_max, st_min = int(words[3]), int(words[1])
                if st_min < n_bp < st_max:
                    print("ERROR in strand divisions!")
                    break
                pad = 3 if st_min <= n_bp else 4
                loc_strand_resid = [i for i in range(st_min, st_max + 1)]
                loc_strand_beads = [i*3-pad+j for i in loc_strand_resid for j in range(3)]
                loc_strand_beads.pop(0)
                strand_resids.append(loc_strand_resid[:])
                strand_beads.append(loc_strand_beads[:])

    gen_pdb_file()
    gen_ninfo_file()
