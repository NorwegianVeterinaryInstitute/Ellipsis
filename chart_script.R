library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

grViz("digraph {
      graph [layout = dot]
      
      node [shape = rectangle, style = filled, fillcolor = '#F6E2BC']
      
      assembly [label = 'Assembly', shape = diamond, fillcolor = '#75BDE0']
      ref [label = 'Plasmid\nreference', shape = diamond, fillcolor = '#75BDE0']
      hybrid [label = 'Hybrid Assembly', shape = diamond, fillcolor = '#75BDE0']
      reads [label = 'Short reads', shape = octagon, fillcolor = '#75BDE0']
      longreads [label = 'Long reads', shape = octagon, fillcolor = '#75BDE0']
      plasmid [label = 'Plasmid fasta', shape = ellipse, fillcolor = '#B9CC95']
      chrom [label = 'Chromosome fasta', shape = ellipse, fillcolor = '#B9CC95']

      fastqc [label = 'FastQC']
      longqc [label = 'Longread QC']
      mobsuite [label = 'Mob-Recon']
      resfinder [label = 'Resfinder']
      virfinder [label = 'Virulencefinder']
      plasfinder [label = 'Plasmidfinder']
      prokka [label = 'Prokka']
      trim [label = 'Trim-galore']
      unicycler [label = 'Unicycler']
      ariba [label = 'ARIBA']
      report [label = 'Collate reports']
      quast [label = 'QUAST']
      filtlong [label = 'Filtlong']
      canu [label = 'Canu']
      bwa [label = 'BWA mem']
      sam [label = 'Samtools']
      cov [label = 'Coverage calc']
      collate [label = 'Collate report']
      
      plasreport [label = 'Plasmid Coverage\nReport', shape = folder, fillcolor = '#E984A2']
      tot [label = 'Total report', shape = folder, fillcolor = '#E984A2']
      res [label = 'Resistance report', shape = folder, fillcolor = '#E984A2']
      vir [label = 'Virulence report', shape = folder, fillcolor = '#E984A2']
      plas [label = 'Plasmid-inc report', shape = folder, fillcolor = '#E984A2']
      sum [label = 'Summary report', shape = folder, fillcolor = '#E984A2']
      
      mobsuite -> plasmid -> resfinder -> report
      longreads -> canu -> filtlong -> unicycler -> hybrid -> mobsuite [color = '#EE8980']
      longreads -> longqc [color = '#EE8980']
      plasmid -> virfinder -> report
      plasmid -> plasfinder -> report
      plasmid -> prokka -> report
      mobsuite -> chrom -> resfinder
      chrom -> virfinder
      chrom -> prokka
      chrom -> plasfinder
      reads -> trim -> unicycler -> assembly -> mobsuite [color = '#1f78b4']
      assembly -> quast [color = '#1f78b4']
      hybrid -> quast [color = '#EE8980']
      reads -> trim -> unicycler [color = '#EE8980']
      reads -> ariba [color = '#1f78b4']
      reads -> ariba [color = '#EE8980']
      report -> tot
      report -> res
      report -> vir
      report -> plas
      report -> sum
      quast -> report
      ariba -> report
      reads -> fastqc [color = '#1f78b4']
      reads -> fastqc [color = '#EE8980']
      
      reads -> bwa -> sam -> cov -> collate -> plasreport [color = '#6a3d9a']
      ref -> bwa [color = '#6a3d9a']
}") %>%
  export_svg() %>% charToRaw() %>% rsvg_png("pipeline.png",
                                            width = 1600,
                                            height = 1700)
