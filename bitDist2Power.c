//-------- BitDist2Power
unsigned int    bitDist2Power(
	unsigned int	*bitDist,		// IN: Bit Distribution in 4-bit (16-level)
	float			*bitPower)		// OUT:Polaris Power
{
	float			wt4[] = {-7.5, -6.5, -5.5, -4.5, -3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
	unsigned int	totalSample = 0;
	int				index_level;

	*bitPower = 0.0;
	for(index_level=0; index_level<16; index_level++){
		totalSample += bitDist[index_level];
		*bitPower	+= wt4[index_level]* wt4[index_level]* (float)bitDist[index_level];
	}
	*bitPower /= (float)totalSample;
	return(totalSample);
}
