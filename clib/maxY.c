float maxY(maxy)
float maxy;
{
	float reduce();
	int a;
	float f;
	float red;
	static float index[] = {1.0,2.0,3.0,5.0,7.0,10.0};
	int counter;
	a = log10( (double) maxy);
	if (maxy < 1.0)
	    a -= 1.0;
	f = pow( (double) 10.0, (double) a);
	for (counter = 0; counter < 6; ++counter) {
	    if (index[counter] * f > maxy)
		break;
	}
	return( index[counter] * f);
}

