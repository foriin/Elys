#!/bin/sh

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/genes.elys.npc.bed --scoreFileName ../bed/kc.nup98.npc.damid.300bp.bigWig --outFileName nup98/metagene.elys.npc.nup98.npc.300.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/genes.elys.nuc.bed --scoreFileName ../bed/kc.nup98.npc.damid.300bp.bigWig --outFileName nup98/metagene.elys.nuc.nup98.npc.300.matrix.gz


plotProfile --startLabel start --endLabel end --averageType mean --matrixFile nup98/metagene.elys.npc.nup98.npc.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName nup98/metaplot.elys.npc.nup98.npc.300.pdf --outFileNameData nup98/metaplot.elys.npc.nup98.npc.300.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile nup98/metagene.elys.nuc.nup98.npc.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName nup98/metaplot.elys.nuc.nup98.npc.300.pdf --outFileNameData nup98/metaplot.elys.nuc.nup98.npc.300.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName genes.elys.bed --scoreFileName ../bed/kc.nup98.npc.damid.300bp.bigWig --outFileName nup98/metagene.elys.other.nup98.npc.300.matrix.gz

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile nup98/metagene.elys.other.nup98.npc.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName nup98/metaplot.elys.other.nup98.npc.300.pdf --outFileNameData nup98/metaplot.elys.other.nup98.npc.300.tab


computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/genes.elys.npc.bed --scoreFileName ../bed/kc.nup98.nucl.damid.300bp.bigWig --outFileName nup98/metagene.elys.npc.nup98.nuc.300.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/genes.elys.nuc.bed --scoreFileName ../bed/kc.nup98.nucl.damid.300bp.bigWig --outFileName nup98/metagene.elys.nuc.nup98.nuc.300.matrix.gz


plotProfile --startLabel start --endLabel end --averageType mean --matrixFile nup98/metagene.elys.npc.nup98.nuc.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName nup98/metaplot.elys.npc.nup98.nuc.300.pdf --outFileNameData nup98/metaplot.elys.npc.nup98.nuc.300.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile nup98/metagene.elys.nuc.nup98.nuc.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName nup98/metaplot.elys.nuc.nup98.nuc.300.pdf --outFileNameData nup98/metaplot.elys.nuc.nup98.nuc.300.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName genes.elys.bed --scoreFileName ../bed/kc.nup98.nucl.damid.300bp.bigWig --outFileName nup98/metagene.elys.other.nup98.nuc.300.matrix.gz

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile nup98/metagene.elys.other.nup98.nuc.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName nup98/metaplot.elys.other.nup98.nuc.300.pdf --outFileNameData nup98/metaplot.elys.other.nup98.nuc.300.tab
