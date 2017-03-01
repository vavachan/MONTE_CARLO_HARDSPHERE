int find_cube(int nAtoms);
void inipos(atom Atoms[],int nAtoms,Vector box);
Vector inipos_bcc(atom Atoms[],int nAtoms,Vector box,double a);
void random_ini(atom Atoms[],int nAtoms,Vector box,double);
Vector FCC(atom Atoms[],int nAtoms,Vector box,double a);
