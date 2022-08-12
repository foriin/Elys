#!/bin/sh

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/genes.elys.npc.bed --scoreFileName ../bed/emb.gaf.chip.300bp.bigWig --outFileName metagene.elys.npc.gaf.300.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/genes.elys.nuc.bed --scoreFileName ../bed/emb.gaf.chip.300bp.bigWig --outFileName metagene.elys.nuc.gaf.300.matrix.gz


plotProfile --startLabel start --endLabel end --averageType mean --matrixFile metagene.elys.npc.gaf.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName metaplot.elys.npc.gaf.300.pdf --outFileNameData metaplot.elys.npc.gaf.300.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile metagene.elys.nuc.gaf.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName metaplot.elys.nuc.gaf.300.pdf --outFileNameData metaplot.elys.nuc.gaf.300.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName genes.elys.bed --scoreFileName ../bed/emb.gaf.chip.300bp.bigWig --outFileName metagene.elys.other.gaf.300.matrix.gz

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile metagene.elys.other.gaf.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName metaplot.elys.other.gaf.300.pdf --outFileNameData metaplot.elys.other.gaf.300.tab


