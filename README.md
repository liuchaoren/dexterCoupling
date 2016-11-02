Tutorial of calculating dexter coupling between pi orbitals, DNA is an example

1. QM calculation (for example with B3LYP/6-31G) of the whole molecule, extract overlap matrix and fockmatrix with rwfdump. Use g09Extract/fockmatrix.py to convert them into numpy array format.

2. With overlap.npy and fockmatrix.npy, use g09Extract/Par_diag.py to block diagonalize the system and obtain homoevecs.npy, lumoevecs.npy which are localized homos and localized lumos in atomic basis. We can also obtain the fockmatrix in localized HOMO, LUMO basis.

3. With the homoevecs.npy and lumoevecs.npy, use PyQuante to calculate the two electron integral between the localized HOMO LUMO basis.

4. With overlap.npy (atomic basis), homoevecs.npy and lumoevecs.npy,  obtain the overlap_HOMOLUMO.npy which is the overlap matrix between localized HOMO and LUMO basis.

5. With overlap_HOMOLUMO.npy, use nonorthogonality/Gram-Schmidt.py to orthogonalize the localized HOMO, LUMO basis.

6. Convert fockmatrix and two electron integrals into the basis of orthogonalized localized HOMO, LUMO basis. 

7. Calculate the H_{ix,jy} with PNAS eq. 3.