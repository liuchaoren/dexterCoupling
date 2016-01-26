'''
calculate two E integrals of four contracted gaussian functions
Author - Chaoren Liu
Date - Jan 25, 2016
'''





def coulomb(a, b, c, d, coul_func):
    Jij = coul_func(a.pexps,a.pcoefs,a.pnorms,a.origin,a.powers,
                    b.pexps,b.pcoefs,b.pnorms,b.origin,b.powers,
                    c.pexps,c.pcoefs,c.pnorms,c.origin,c.powers,
                    d.pexps,d.pcoefs,d.pnorms,d.origin,d.powers)
    return a.norm*b.norm*c.norm*d.norm*Jij





