#!/bin/sh

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/silent.genes.elys.x.npc.bed --scoreFileName ../bed/s2.gaf.chip.300bp.bigWig --outFileName silact/metagene.silent.elys.npc.s2.gaf.300.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/active.genes.elys.x.npc.bed --scoreFileName ../bed/s2.gaf.chip.300bp.bigWig --outFileName silact/metagene.active.elys.npc.s2.gaf.300.matrix.gz


plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.silent.elys.npc.s2.gaf.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.silent.elys.npc.s2.gaf.300.pdf --outFileNameData silact/metaplot.silent.elys.npc.s2.gaf.300.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.active.elys.npc.s2.gaf.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.active.elys.npc.s2.gaf.300.pdf --outFileNameData silact/metaplot.active.elys.npc.s2.gaf.300.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/silent.genes.elys.x.npc.bed --scoreFileName ../bed/emb.gaf.chip.300bp.bigWig --outFileName silact/metagene.silent.elys.npc.gaf.300.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/active.genes.elys.x.npc.bed --scoreFileName ../bed/emb.gaf.chip.300bp.bigWig --outFileName silact/metagene.active.elys.npc.gaf.300.matrix.gz


plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.silent.elys.npc.gaf.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.silent.elys.npc.gaf.300.pdf --outFileNameData silact/metaplot.silent.elys.npc.gaf.300.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.active.elys.npc.gaf.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.active.elys.npc.gaf.300.pdf --outFileNameData silact/metaplot.active.elys.npc.gaf.300.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/silent.genes.elys.x.npc.bed --scoreFileName ../bed/s2.pol2.chip.300bp.bigWig --outFileName silact/metagene.silent.elys.npc.pol2.300.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/active.genes.elys.x.npc.bed --scoreFileName ../bed/s2.pol2.chip.300bp.bigWig --outFileName silact/metagene.active.elys.npc.pol2.300.matrix.gz

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.silent.elys.npc.pol2.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.silent.elys.npc.pol2.300.pdf --outFileNameData silact/metaplot.silent.elys.npc.pol2.300.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.active.elys.npc.pol2.300.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.active.elys.npc.pol2.300.pdf --outFileNameData silact/metaplot.active.elys.npc.pol2.300.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/silent.genes.elys.x.npc.bed --scoreFileName dm3.300bp.at.perc.bigWig --outFileName silact/metagene.silent.elys.npc.at.count.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/active.genes.elys.x.npc.bed --scoreFileName dm3.300bp.at.perc.bigWig --outFileName silact/metagene.active.elys.npc.at.count.matrix.gz

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.silent.elys.npc.at.count.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.silent.elys.npc.at.count.pdf --outFileNameData silact/metaplot.silent.elys.npc.at.count.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.active.elys.npc.at.count.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.active.elys.npc.at.count.pdf --outFileNameData silact/metaplot.active.elys.npc.at.count.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/silent.genes.elys.x.npc.bed --scoreFileName elys.embryos.bigWig --outFileName silact/metagene.silent.elys.npc.emb.elys.pr.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/active.genes.elys.x.npc.bed --scoreFileName elys.embryos.bigWig --outFileName silact/metagene.active.elys.npc.emb.elys.pr.matrix.gz

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.silent.elys.npc.emb.elys.pr.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.silent.elys.npc.emb.elys.pr.pdf --outFileNameData silact/metaplot.silent.elys.npc.emb.elys.pr.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.active.elys.npc.emb.elys.pr.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.active.elys.npc.emb.elys.pr.pdf --outFileNameData silact/metaplot.active.elys.npc.emb.elys.pr.tab

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/silent.genes.elys.x.npc.bed --scoreFileName ../bed/s2.h3k27ac.chip.300bp.bigWig --outFileName silact/metagene.silent.elys.npc.s2.h3k27ac.matrix.gz

computeMatrix scale-regions -p 10 --startLabel start --endLabel end --upstream 500 --downstream 500 --regionBodyLength 5000 --regionsFileName ../bed/active.genes.elys.x.npc.bed --scoreFileName ../bed/s2.h3k27ac.chip.300bp.bigWig --outFileName silact/metagene.active.elys.npc.s2.h3k27ac.matrix.gz

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.silent.elys.npc.s2.h3k27ac.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.silent.elys.npc.s2.h3k27ac.pdf --outFileNameData silact/metaplot.silent.elys.npc.s2.h3k27ac.tab

plotProfile --startLabel start --endLabel end --averageType mean --matrixFile silact/metagene.active.elys.npc.s2.h3k27ac.matrix.gz --plotHeight 15 --plotWidth 20 --outFileName silact/metaplot.active.elys.npc.s2.h3k27ac.pdf --outFileNameData silact/metaplot.active.elys.npc.s2.h3k27ac.tab
