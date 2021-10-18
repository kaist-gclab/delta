/*
qsmodel.h     headerfile for quasistatic probability model

(c) Michael Schindler
1997, 1998, 2000
http://www.compressconsult.com/
michael@compressconsult.com

c++ port by Sebastien Valette

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.  It may be that this
program violates local patents in your country, however it is
belived (NO WARRANTY!) to be patent-free here in Austria.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA.

Qsmodel is a quasistatic probability model that periodically
(at chooseable intervals) updates probabilities of symbols;
it also allows to initialize probabilities. Updating is done more
frequent in the beginning, so it adapts very fast even without
initialisation.

it provides function for creation, deletion, query for probabilities
and symbols and model updating.

*/

#ifndef QUASISTATICMODEL_H
#define QUASISTATICMODEL_H
#include <iostream>

#include "port.h"
#include <vtkCommand.h>

class VTK_EXPORT qsmodel 
{
public:

	/* initialisation of qsmodel                           */
	/* m   qsmodel to be initialized                       */
	/* n   number of symbols in that model                 */
	/* lg_totf  base2 log of total frequency count         */
	/* rescale  desired rescaling interval, should be < 1<<(lg_totf+1) */
	/* init  array of int's to be used for initialisation (NULL ok) */
	/* compress  set to 1 on compression, 0 on decompression */
	void initqsmodel(int n, int lg_totf, int rescale,
		int *init, int compress );

	/* reinitialisation of qsmodel                         */
	/* m   qsmodel to be initialized                       */
	/* init  array of int's to be used for initialisation (NULL ok) */
	void resetqsmodel(int *init);

	/* retrieval of estimated frequencies for a symbol     */
	/* m   qsmodel to be questioned                        */
	/* sym  symbol for which data is desired; must be <n   */
	/* sy_f frequency of that symbol                       */
	/* lt_f frequency of all smaller symbols together      */
	/* the total frequency is 1<<lg_totf                   */
	void qsgetfreq(int sym, int *sy_f, int *lt_f );


	/* find out symbol for a given cumulative frequency    */
	/* m   qsmodel to be questioned                        */
	/* lt_f  cumulative frequency                          */
	int qsgetsym(int lt_f );

	/* update model                                        */
	/* m   qsmodel to be updated                           */
	/* sym  symbol that occurred (must be <n from init)    */
	void qsupdate(int sym );

	qsmodel();
	~qsmodel();

private:

	/* rescale frequency counts */
	void dorescale();

	/* deletion of qsmodel m                               */
	void deleteqsmodel();

	int n;           /* number of symbols */
	int left;          /* symbols to next rescale */
	int nextleft;      /* symbols with other increment */
	int rescale;       /* intervals between rescales */
	int targetrescale; /* should be interval between rescales */
	int incr;          /* increment per update */
	int searchshift;   /* shift for lt_freq before using as index */
	uint4 *cf;         /* array of cumulative frequencies */
	int *newf;         /* array for collecting ststistics */
	int *search;       /* structure for searching on decompression */
};
#endif
