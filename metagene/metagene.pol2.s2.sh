#!/bin/sh

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/genes.elys.npc.bed --scoreFileName ../bed/s2.pol2.chip.300bp.bigWig --outFileName metagene.elys.npc.pol2.300.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/genes.elys.nuc.bed --scoreFileName ../bed/s2.pol2.chip.300bp.bigWig --outFileName metagene.elys.nuc.pol2.300.matrix.gz


plotProfile --startLabel start --endLabel end --averageType mean --matrixFile metagene.elys.npc.pol2.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName metaplot.elys.npc.pol2.300.pdf --outFileNameData metaplot.elys.npc.pol2.300.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile metagene.elys.nuc.pol2.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName metaplot.elys.nuc.pol2.300.pdf --outFileNameData metaplot.elys.nuc.pol2.300.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName genes.elys.bed --scoreFileName ../bed/s2.pol2.chip.300bp.bigWig --outFileName metagene.elys.other.pol2.300.matrix.gz

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile metagene.elys.other.pol2.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName metaplot.elys.other.pol2.300.pdf --outFileNameData metaplot.elys.other.pol2.300.tab


