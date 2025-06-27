### Prototype

The prototype of the REMAP pipeline is desgined to have multiple rounds of mappings.

In each round, the pipeline detects the 5' softclipped reads, and remove the last unmapped nucleotide.

The trimmed reads were then remapped to the reference. 

The mapping cycles will be performed for 10 times.

However, I found that this method cannot precisely reveal expansion events. Based on the overall understanding of RNA 5' epansions, I designed the current pipeline.