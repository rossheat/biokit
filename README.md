<div style="text-align:center"><img src="https://biokit.org/img/banner.png" /></div>

BioKit is a Dart package for Bioinformatics.

Ensure that you have BioKit [installed](https://biokit.org/docs/) before continuing.

This document is intended to make you proficient with BioKit in the least amount of time possible; you can read through it sequentially, or if you're reading this on [biokit.org](https://biokit.org/docs/usage), use the heading menu on the right side of the page to jump to a topic of interest. 

If you want a deeper look at how BioKit works, view our [API Reference](https://pub.dev/documentation/biokit/latest/biokit/biokit-library.html).

## Creating Sequences

Create a `DNA`, `RNA` or `Peptide` instance: 

```dart
DNA dnaSeq = DNA(seq: 'ATGCTA');

RNA rnaSeq = RNA(seq: 'AUGCUA');

Peptide pepSeq = Peptide(seq: 'MSLAKR');
```

`DNA` and `RNA` classes must be initialized with a `String` of at least six valid nucleotides, while the `Peptide` class requires a minimum of two valid amino acids. 

If any monomer in the sequence passed to the `seq` parameter is not valid for the class, an `error` is thrown. 

### Add Sequence Metadata

Optionally, you can add `name`, `id`, and `desc` metadata when you instantiate the class. Using `DNA` as an example: 

```dart
DNA dnaSeq = DNA(seq: 'ATGCTA', name: 'My Name', id: 'My ID', desc: 'My Description');
```

If you do not set a value for the `name`, `id`, or `desc` fields at the time of instantiation, each will receive a default `String` value. 

## Get Properties

Return the values of the properties of a `DNA`, `RNA`, or `Peptide` instance: 

```dart
dnaSeq.seq;
// ATGCTA

dnaSeq.len;
// 6 

dnaSeq.id;
// Default ID

dnaSeq.name;
// Default name

dnaSeq.desc;
// Default description

dnaSeq.type;
// dna
```

## Set Properties

Update the properties of a `DNA`, `RNA`, or `Peptide` instance:

```dart
dnaSeq.name = 'New name';

dnaSeq.id = 'New ID';

dnaSeq.desc = 'New description';
```

## Sequence Info

View information about a `DNA`, `RNA`, or `Peptide` instance by calling its `info()` method or printing it to the console:

```dart
dnaSeq.info();
/*
{
   "seq":"ATGCTA",
   "type":"dna",
   "monomers":6,
   "name":"New Name",
   "id":"New ID",
   "desc":"New description"
}
*/

print(dnaSeq);
/*
{
   "seq":"ATGCTA",
   "type":"dna",
   "monomers":6,
   "name":"New Name",
   "id":"New ID",
   "desc":"New description"
}
*/
```

## Random Sequences

Return a random `DNA`, `RNA`, or `Peptide` instance with the `random()` method and pass the desired length of the sequence to the `len` parameter: 

```dart
// A random DNA instance with 20 nucleotides.
DNA dnaSeq = DNA.random(len: 20);

dnaSeq.info();

/*
{
   "seq":"TAACTTCGATCGCTCTGGCA",
   "type":"dna",
   "monomers":20,
   "name":"Default Name",
   "id":"Default ID",
   "desc":"Default description"
}
*/
```

## FASTA Data

BioKit contains a number of methods and functions for working with FASTA formatted data.

### Uniprot ID

Return a `String` of protein data in FASTA format using the static `uniprotIdToFASTA()` method from the `Utils` class:

```dart

String proteinFASTA = await Utils.uniprotIdToFASTA(uniprotId: 'B5ZC00');

/*
>sp|B5ZC00|SYG_UREU1 Glycine--tRNA ligase OS=Ureaplasma urealyticum ...
MKNKFKTQEELVNHLKTVGFVFANSEIYNGLANAWDYGPLGVLLKNNLKNLWWKEFVTKQ
KDVVGLDSAIILNPLVWKASGHLDNFS ...
*/
```

Note that this method requires network access. 

### Read `String`

Use the `readFASTA()` method to parse FASTA formatted `String` data. 

`readFASTA()` is able to parse FASTA files containing multiple sequences, and hence returns a `List`:

```dart
List<Map<String, String>> proteinMaps = await Utils.readFASTA(str: proteinFASTA);

/*
[
   {
      "seq":"MKNKFKTQEELVNHLKTVGFVFANSEIYNGLANAWDYGPLGVLLKNNLKNLWWKEFVTK ... ",
      "id":"sp|B5ZC00|SYG_UREU1",
      "desc":"Glycine--tRNA ligase OS=Ureaplasma urealyticum serovar 10 (... "
   }
]
*/
```

### Read File

Read in data from a FASTA formatted txt file:

```dart
List<Map<String, String>> dnaMaps = await Utils.readFASTA(path: './gene_bank.txt');

/*
[
   {
      "seq":"GGCAGATTCCCCCTAGACCCGCCCGCACCATGGTCAGGCATGCCCCTCCTCATCGCTGG ... ",
      "id":"HSBGPG",
      "desc":"Human gene for bone gla protein (BGP)"
   },
   {
      "seq":"CCACTGCACTCACCGCACCCGGCCAATTTTTGTGTTTTTAGTAGAGACTAAATACCATA ... ",
      "id":"HSGLTH1",
      "desc":"Human theta 1-globin gene"
   }
]
*/
```

### Write File

Write the contents of a `DNA`, `RNA`, or `Peptide` instance to a FASTA formatted txt file using the `toFASTA()` method:

```dart
// Get the first Map object.
Map<String, String> firstSeq = dnaMaps.first;

// Create a new DNA instance.
DNA dnaSeq = DNA(seq: firstSeq['seq']!, id: firstSeq['id']!, desc: firstSeq['desc']!);

// Write the instance contents to FASTA formatted file.
dnaSeq.toFASTA(path: '../deliverables', filename: 'my_dna_seq');

/*
>HSBGPG Human gene for bone gla protein (BGP)
GGCAGATTCCCCCTAGACCCGCCCGCACCATGGTCAGGCATGCCCCTCCTCATCGCTGGG
CACAGCCCAGAGGGTATAAACAGTGCTGGAGGCTGGCGGGGCAGGCCAGCTGAGTCCTGA
GCAGCAGCCCAGCGCAGCCACCGAGACA ...
*/
```

## DNA Analysis Report

Create a DNA analysis report by calling the `report()` method on a `DNA` instance:

```dart
dnaSeq.report(path: '../deliverables', creator: 'John Doe', title: 'BGP Report');
```

## + Operator

Return the concatenated sequence result of two or more `DNA`, `RNA`, or `Peptide` instance sequences, of the same type, with the `+` operator:

```dart
RNA rnaSeq1 = RNA(seq: 'AUGCAG');
RNA rnaSeq2 = RNA(seq: 'GCUGAA');

rnaSeq1 + rnaSeq2; 
// "AUGCAGGCUGAA"
```

## Reversing

Reverse a `DNA`, `RNA`, or `Peptide` instance's sequence with the `reverse()` method: 

```dart
Peptide pepSeq = Peptide(seq: 'MPAG');

pepSeq.reverse();
// GAPM
```

## Point Mutations

Return the number of positional-differences between two `DNA`, `RNA`, or `Peptide`  instance sequences, of the same type, with the `difference()` method:

```dart
DNA dnaSeq1 = DNA(seq: 'ATGCAT');

// Difference: "A" at index 1, and "T" at index 4. 
DNA dnaSeq2 = DNA(seq: 'AAGCTT');

dnaSeq1.difference(oSeq: dnaSeq2)
// 2
```

## Motif Detection

BioKit has a number of functions and methods to convert and detect matches between a motif and the sequence of a `DNA`, `RNA`, or `Peptide` instance. 

### Find Motifs

Return the indices of all matches between a `DNA`, `RNA`, or `Peptide` instance's sequence and the sequence passed to the `findMotif()` method's `motif` parameter:

```dart
RNA rnaSeq = RNA(seq: 'GAUAUAUC');

rnaSeq.findMotif(motif: 'AUAU');

/*
{
   "matchCount":2,
   "matchIndices":[
      {
         "match":"AUAU",
         "startIndex":1,
         "endIndex":4
      },
      {
         "match":"AUAU",
         "startIndex":3,
         "endIndex":6
      }
   ]
}
*/
```

Set `overlap` to `false` to return only the match indices that do not overlap: 

```dart
rnaSeq.findMotif(motif: 'AUAU', overlap: false);

/*
{
   "matchCount":1,
   "matchIndices":[
      {
         "match":"AUAU",
         "startIndex":0,
         "endIndex":3
      }
   ]
}
*/
```

### Shared Motifs

Return the longest shared motif between two `DNA`, `RNA`, or `Peptide` instance sequences, of the same type: 

```dart
DNA dnaSeq1 = DNA('GATATA');

DNA dnaSeq2 = DNA('AGCATA');

dnaSeq1.sharedMotif(oSeq: dnaSeq2); 
// ATA
```

### Manually Convert Motif to Regex

The `findMotif()` method automatically converts motifs passed to its `motif` parameter to regular-expression format, however, you can also perform the conversion manually using the `motifToRe()` function: 

```dart
Utils.motifToRe(motif: 'N{P}[ST]{P}'); 
// 'N[^P][S|T|][^P]'

// No change needs to be made.
Utils.motifToRe(motif: 'ATGC');
// ATGC
```

## Splicing

Return a sequence with all occurrences of a motif removed from a `DNA`, `RNA`, or `Peptide` instance's sequence using the `splice` method, and passing the motif to the `motif` parameter: 

```dart
RNA rnaSeq = RNA(seq: 'AUCAUGU');

// Removes all occurrences of 'AU'.
rnaSeq.splice(motif: 'AU');
// CGU
```

## Monomer Frequency

Return the frequency of each monomer in a `DNA`, `RNA`, or `Peptide` instance's sequence with the `freq()` method:

```dart
DNA dnaSeq = DNA(seq: 'AGCTTTTCAGC');

dnaSeq.freq();

/*
{
   "A":2.0,
   "G":2.0,
   "C":3.0,
   "T":4.0
}
*/
```

### Percentage of Total

Return the percentage of the total that each monomer count represents in the sequence by passing `true` to the `norm` parameter of the `freq()` method:

```dart
dnaSeq.freq(norm: true);

/*
{
   "A":18.2,
   "G":18.2,
   "C":27.3,
   "T":36.4
}
*/
```

### Ignore the Stop Amino Acid

When the `translate()` method is called on `DNA` or `RNA` instances, BioKit returns an amino acid sequence; when BioKit encounters a stop codon, rather than stoping translation, or ignoring the stop codon, BioKit places an "X" character at that position in the amino acid sequence:

```dart
// UAG is a stop codon
RNA rnaSeq = RNA(seq: 'CGGUAGACU'); 

rnaSeq.translate();

/*
{
   "aaSeq":"RXT",
   "nucCount":8,
   "aaCount":3
}
*/
```

Therefore, If you use the `aaSeq` key's value to create a new `Peptide` instance, and then execute the `freq()` method, the "X" character will be taken into account as part of the calculation:

```dart
// Create a Peptide instance using the RNA instance translation product.
Peptide pepSeq = Peptide(seq: rnaSeq.translate()['aaSeq']!);

pepSeq.freq(); 

/*
{
   "R":1.0,
   "X":1.0,
   "T":1.0
}
*/ 
```

However, if you do not want the "X" character to be taken into account as part of the calculation, pass `true` to the `ignoreStopAA` parameter of the `freq()` method:

```dart
pepSeq.freq(ignoreStopAA: true);

/*
{
   "R":1.0,
   "T":1.0
}
*/
```

## Modified Sequence Length

In addition to being able to return the length of a  `DNA`, `RNA`, or `Peptide` instance's sequence by using the `len` getter:

```dart
DNA dnaSeq = DNA(seq: 'ATGCGAT');

dnaSeq.len;
// 7 
```

You can also return the length of the sequence minus a particular monomer by using the `lenMinus()` method, and passing the `monomer` you'd like to discount:

```dart
dnaSeq.lenMinus(monomer: 'A');
// 5
```

## Generate Combinations

Return all possible combinations of a  `DNA`, `RNA`, or `Peptide` instance's sequence using the `combinations()` method:

```dart
Peptide pepSeq = Peptide(seq: 'MSTC');

pepSeq.combinations(); 
// [M, MS, MST, MSTC, S, ST, STC, T, TC]
```

Sort the combinations by setting `sorted` to `true`:

```dart
pepSeq.combinations(sorted: true);
// [MSTC, MST, STC, MS, ST, TC, M, S, T]
```

## Codon Frequency

Return the frequency of a codon in a `DNA` or `RNA` instance's sequence using the `codonFreq()` method, passing the codon of interest to the `codon` parameter:

```dart
RNA rnaSeq = RNA(seq: 'AUGAGGAUGCACAUG');

rnaSeq.codonFreq(codon: 'AUG');
// 3 
```

Be aware that `codonFreq()` scans the sequence in batches of three nucleotides per step, starting with the first three nucleotides in the sequence. Therefore, the exact `codon` must be present in a batch in order to be detected.

## Complementary Strand

Return the complementary strand to a `DNA` or `RNA` instance sequence's with the `complementary()` method:

```dart
DNA dnaSeq = DNA(seq: 'AAACCCGGT');

dnaSeq.complementary();
// TTTGGGCCA
```

To return the reverse complementary strand, pass `true` to the `rev` parameter:

```dart
dnaSeq.complementary(rev: true);
// ACCGGGTTT
```

## Guanine & Cytosine Content

Return the percentage of Guanine and Cytosine content in a `DNA` or `RNA` instance's sequence with the `gcContent()` method:

```dart
DNA dnaSeq = DNA(seq: 'TCCCTACGCCG');

dnaSeq.gcContent();
// 72.73
```

## Translation

Return the amino acid translation product from a `DNA` or `RNA` instance's sequence, using the `translate()` method: 

```dart
RNA rnaSeq = RNA(seq: 'AUGGCCAUGGCGCCCAGAACU');

rnaSeq.translate();

/*
{
   "aaSeq":"MAMAPRT",
   "nucCount":20,
   "aaCount":7
}
*/
```

Return the reverse complementary translation strand by passing `true` to the `rev` parameter: 

```dart
rnaSeq.translate(rev: true); 

/*
{
   "aaSeq":"SSGRHGH",
   "nucCount":20,
   "aaCount":7
}
*/
```

Modify the index in which translation starts by passing the desired start index to the `startIdx` parameter:

```dart
rnaSeq.translate(startIdx: 2);

/*
{
   "aaSeq":"GHGAQN",
   "nucCount":18,
   "aaCount":6
}
*/
```

## Generate Proteins

Return proteins from open reading frames present in a `DNA` or `RNA` instance sequence's with the `proteins()` method: 

```dart
DNA dnaSeq = DNA(seq: 'AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCTGAATGATCCGAGTAGCATCTCAG');

dnaSeq.proteins(); 
// [MLLGSFRGHPHVT, MGMTPE, MTPE, M, M]
```

Return only unique proteins by passing `true` to the `unique` parameter:

```dart
dnaSeq.proteins(unique: true);
// [MLLGSFRGHPHVT, MGMTPE, MTPE, M]
```

## Transcription

Return the RNA transcription product from a `DNA` instance's sequence using the `transcribe()` method:

```dart
DNA dnaSeq = DNA(seq: 'TACGTAA');

dnaSeq.transcribe();
// UACGUAA
```

Change where transcription starts from by passing the desired start index to the `startIdx` parameter:

```dart
dnaSeq.transcribe(startIdx: 3); 
// GUAA
```

## Restriction Sites

Return restriction sites in a `DNA` instance's sequence with the `restrictionSites()` method: 

```dart
DNA dnaSeq = DNA(seq: 'TGCATGTCTATATG');

dnaSeq.restrictionSites();

/*
{
   "TGCA":[
      {
         "startIdx":0,
         "endIndex":4
      }
   ],
   "CATG":[
      {
         "startIdx":2,
         "endIndex":6
      }
   ],
   "TATA":[
      {
         "startIdx":8,
         "endIndex":12
      }
   ],
   "ATAT":[
      {
         "startIdx":9,
         "endIndex":13
      }
   ]
}
*/
```

Pass values to the `minSiteLen` and `maxSiteLen` parameters to change the restriction site search length. 

## Transition/Transversion Ratio

Return the transition/transversion ratio between two `DNA` instance sequences with the `tranRatio()` method:

```dart
DNA dnaSeq1 = DNA(seq: 'GACTGGTGGAAGT');

DNA dnaSeq2 = DNA(seq: 'TTATCGGCTGAAT');

dnaSeq1.tranRatio(oSeq: dnaSeq2); 
// 0.29
```

Note that if the number of transversions is equal to `0`, the method returns `-1`, as division by `0` is undefined and leads to a result of `inf`. 

## Double Helix Geometric Length

Return the geometric length (nm) of a double helix formed by a `DNA` instance's sequence using the `dHelixGeoLen()` method:

```dart
DNA dnaSeq = DNA(seq: 'ATGCATGC');

dnaSeq.dHelixGeoLen();
// 2.72
```

## Double Helix Turns

Return the number of turns in a double helix formed by a `DNA` instance's sequence using the `dHelixTurns()` method:

```dart
DNA dnaSeq = DNA(seq: 'ATGCATGCATGCATGC');

dnaSeq.dHelixTurns();
// 1.6 
```

## Reverse Transcription

Return the reverse transcription product from an `RNA` instance's sequence using the `revTranscribe()` method:

```dart
RNA rnaSeq = RNA(seq: 'AUGCUAGU');

rnaSeq.revTranscribe();
// ATGCTAGT
```

## Monoisotopic Mass

Return the Monoisotopic mass (Da) of a `Peptide` instance's sequence using the `monoMass()` method: 

```dart
Peptide pepSeq = Peptide(seq: 'MSTGARVD');

pepSeq.monoMass();
// 817.38
```

Modify the number of decimal places by passing a the desired number of decimals to the `decimals` parameter:

```dart
pepSeq.monoMass(decimals: 1);
// 817.4
```

Return the Monoisotopic mass in kDa by passing `true` to the `kDa` parameter:

```dart
pepSeq.monoMass(kDa: true);
// 0.82
```

