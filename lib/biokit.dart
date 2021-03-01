// import 'functions/DNAtoRNA.dart';
// import 'functions/GCContent.dart';
// import 'functions/RNAtoProtein.dart';
// import 'functions/complementaryStrand.dart';
// import 'functions/dominantPhenotypeProbability.dart';
// import 'functions/hammingDistance.dart';
// import 'functions/matchMotif.dart';
// import 'functions/nucleotideFrequency.dart';
// import 'helpers/detectSequenceType.dart';
// import 'helpers/makeNucleotideSequence.dart';
// import 'helpers/validateNucleotideSequence.dart';

export 'functions/DNAtoRNA.dart';
export 'functions/GCContent.dart';
export 'functions/RNAtoProtein.dart';
export 'functions/complementaryStrand.dart';
export 'functions/dominantPhenotypeProbability.dart';
export 'functions/hammingDistance.dart';
export 'functions/matchMotif.dart';
export 'functions/nucleotideFrequency.dart';
export 'helpers/detectSequenceType.dart';
export 'helpers/makeNucleotideSequence.dart';
export 'helpers/validateNucleotideSequence.dart';

// void main() {
//   String DNASequence = makeNucleotideSequence(40, 'DNA');
//   String RNASequence = makeNucleotideSequence(40, 'RNA');
//   print(validateNucleotideSequence(DNASequence, 'DNA'));
//   print(nucleotideFrequency(DNASequence, 'DNA'));
//   print(DNAToRNA('GATGGAACTTGACTACGTAAATT'));
//   print(complementaryStrand('AAAACCCGGT', 'DNA', true));
//   print(GCContent('CCACCCTCGTGGTAGGCAGTAGGTGGAAT', 'DNA'));
//   print(hammingDistance('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT', 'DNA'));
//   print(dominantPhenotypeProba(19, 29, 28)); // 0.6892982456140351
//   print(detectSequenceType('AGCGG'));
//   print(RNAtoProtein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'));
//   print(matchMotif('GATATATGCATATACTT', 'ATAT', 'DNA'));
// }
