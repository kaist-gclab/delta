/*
rangecod.c     range encoding

(c) Michael Schindler
1997, 1998, 1999, 2000
http://www.compressconsult.com/
michael@compressconsult.com

 C++ port by Sebastien Valette


This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.  It may be that this
program violates local patents in your country, however it is
belived (NO WARRANTY!) to be patent-free. Glen Langdon also
confirmed my poinion that IBM UK did not protect that method.


You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA.

Range encoding is based on an article by G.N.N. Martin, submitted
March 1979 and presented on the Video & Data Recording Conference,
Southampton, July 24-27, 1979. If anyone can name the original
copyright holder of that article please contact me; this might
allow me to make that article available on the net for general
public.

Range coding is closely related to arithmetic coding, except that
it does renormalisation in larger units than bits and is thus
faster. An earlier version of this code was distributed as byte
oriented arithmetic coding, but then I had no knowledge of Martin's
paper from 1997.

The input and output is done by the inbyte and outbyte macros
defined in the .c file; change them as needed; the first parameter
passed to them is a pointer to the rangedecoder structure; extend that
structure as needed (and don't forget to initialize the values in
start_encoding resp. start_decoding). This distribution writes to
stdout and reads from stdin.

There are no global or static var's, so if the IO is thread save the
whole rangedecoder is - unless GLOBALRANGECODER in rangecod.h is defined.

For error recovery the last 3 bytes written contain the total number
of bytes written since starting the encoder. This can be used to
locate the beginning of a block if you have only the end. The header
size you pass to initrangecoder is included in that count.

There is a supplementary file called renorm95.c available at the
website (www.compressconsult.com/rangecoder/) that changes the range
coder to an arithmetic coder for speed comparisons.*/

/*
define NOWARN if you do not expect more than 2^32 outstanding bytes 
since I recommend restarting the coder in intervals of less than    
2^23 symbols for error tolerance this is not expected
*/
#define NOWARN

/*
define EXTRAFAST for increased speed; you loose compression and
compatibility in exchange.
*/
/* #define EXTRAFAST */

#include "port.h"
#include "RangeDecoder.h"

/* SIZE OF RANGE ENCODING CODE VALUES. */

#define CODE_BITS 32
#define Top_value ((code_value)1 << (CODE_BITS-1))

#define SHIFT_BITS (CODE_BITS - 9)
#define EXTRA_BITS ((CODE_BITS-2) % 8 + 1)
#define Bottom_value (Top_value >> 8)

unsigned char rangedecoder::inbyte()
{
	unsigned char TempChar;
	this->InputFile.read((char*) &TempChar,sizeof (char));
	return (TempChar);
}


void rangedecoder::open_file(const char* file)
{ 
	this->InputFile.open (file, std::ofstream::in|std::ios::binary);
}
void rangedecoder::close_file()
{ 
	this->InputFile.close();
};

/* Start the decoder                                         */
/* rc is the range coder to be used                          */
/* returns the char from start_encoding or EOF               */
int rangedecoder::start_decoding()
{
	char c=this->inbyte();
	if (c==EOF)
		return EOF;
	this->buffer = this->inbyte();
	this->low = this->buffer >> (8-EXTRA_BITS);
	this->range = (code_value)1 << EXTRA_BITS;
	return c;
}


void rangedecoder::dec_normalize()
{  
	while (this->range <= Bottom_value)
	{   
		this->low = (this->low<<8) | ((this->buffer<<EXTRA_BITS)&0xff);
		this->buffer = this->inbyte();
		this->low |= this->buffer >> (8-EXTRA_BITS);
		this->range <<= 8;
	}
}

/* Calculate culmulative frequency for next symbol. Does NO update!*/
/* rc is the range coder to be used                          */
/* tot_f is the total frequency                              */
/* or: totf is (code_value)1<<shift                                      */
/* returns the culmulative frequency                         */
freq rangedecoder::decode_culfreq(freq tot_f )
{  
	freq tmp;
	this->dec_normalize();
	this->help = this->range/tot_f;
	tmp = this->low/this->help;
#ifdef EXTRAFAST
	return tmp;
#else
	return (tmp>=tot_f ? tot_f-1 : tmp);
#endif
}

freq rangedecoder::decode_culshift(freq shift )
{ 
	freq tmp;
	this->dec_normalize();
	this->help = this->range>>shift;
	tmp = this->low/this->help;
#ifdef EXTRAFAST
	return tmp;
#else
	return (tmp>>shift ? ((code_value)1<<shift)-1 : tmp);
#endif
}


/* Update decoding state                                     */
/* rc is the range coder to be used                          */
/* sy_f is the interval length (frequency of the symbol)     */
/* lt_f is the lower end (frequency sum of < symbols)        */
/* tot_f is the total interval length (total frequency sum)  */
void rangedecoder::decode_update(freq sy_f, freq lt_f, freq tot_f)
{   
	code_value tmp;
	tmp = this->help * lt_f;
	this->low -= tmp;
#ifdef EXTRAFAST
	this->range = this->help * sy_f;
#else
	if (lt_f + sy_f < tot_f)
		this->range = this->help * sy_f;
	else
		this->range -= tmp;
#endif
}


/* Decode a byte/short without modelling                     */
/* rc is the range coder to be used                          */
unsigned char rangedecoder::decode_byte()
{   
	unsigned char tmp = this->decode_culshift(8);
	this->decode_update(1,tmp,(freq)1<<8);
	return tmp;
}

unsigned short rangedecoder::decode_short()
{  
	unsigned short tmp = this->decode_culshift(16);
	this->decode_update(1,tmp,(freq)1<<16);
	return tmp;
}


/* Finish decoding                                           */
/* rc is the range coder to be used                          */
void rangedecoder::done_decoding()
{  
	this->dec_normalize();      /* normalize to use up all bytes */
}
