int	pow2round(
	unsigned int	number)		// Input Number
{
	unsigned int	temp;
	int				bitCount = 0;
//-------- Loop to Number of bits
	temp = number;
	while(temp > 0){
		temp = temp >> 1;
		bitCount ++;
	}
//-------- Round number
	temp = 0x01 << (bitCount - 1);
	if(temp != number){temp = temp << 1;	}
	return(temp);
}
