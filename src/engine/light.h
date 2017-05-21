#ifndef LIGHT 
#define LIGHT

class light{
public:
	GLenum l;
	float pos[4], amb[4], diff[4], spec[4];
	float dir[3];
	float exp, cut;

	light(GLenum el, float fpos[], float famb[], float fdiff[], float fspec[], float fdir[], float fexp, float fcut) : 
		l{el}, pos{fpos[0], fpos[1], fpos[2], fpos[3]}, amb{famb[0], famb[1], famb[2], famb[3]},
		 diff{fdiff[0], fdiff[1], fdiff[2], fdiff[3]}, spec{fspec[0], fspec[1], fspec[2], fspec[3]},
		 dir{fdir[0], fdir[1], fdir[2]}, exp{fexp}, cut{fcut} {};

	void apply(){
		glLightfv(l, GL_POSITION, pos);
		glLightfv(l, GL_AMBIENT, amb);
		glLightfv(l, GL_DIFFUSE, diff);
		glLightfv(l, GL_SPECULAR, spec);
		glLightfv(l, GL_SPOT_DIRECTION, dir);
		glLightf(l, GL_SPOT_EXPONENT, exp);
		glLightf(l, GL_SPOT_CUTOFF, cut);
	}
};

#endif