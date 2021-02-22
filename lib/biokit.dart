library biokit;

import 'dart:math';

const List<String> DNANucleotides = ["A", "T", "G", "C"];
const List<String> RNANucleotides = ["A", "U", "G", "C"];

const Map<String, String> DNAComplementaryNucleotides = {
  'A': 'T',
  'T': 'A',
  'G': 'C',
  'C': 'G'
};
const Map<String, String> RNAComplementaryNucleotides = {
  'A': 'U',
  'U': 'A',
  'G': 'C',
  'C': 'G'
};

const String kDNA = 'dna';
const String kRNA = 'rna';

const String kSequence = 'sequence';
const String kType = 'type';

Map<String, String> validateNucleotideSequence(
    String nucleotideSequence, String sequenceType) {
  String upperNucleotideSequence = nucleotideSequence.toUpperCase();

  String lowerSequenceType = _validateSequenceType(sequenceType);

  upperNucleotideSequence.split('').asMap().forEach((index, nucleotide) {
    if (lowerSequenceType == kDNA
        ? !DNANucleotides.contains(nucleotide)
        : !RNANucleotides.contains(nucleotide)) {
      throw (_invalidSequenceErrorMessage(
          nucleotide = nucleotide, index = index, sequenceType = sequenceType));
    }
  });

  return {kSequence: upperNucleotideSequence, kType: lowerSequenceType};
}

String _invalidSequenceTypeErrorMessage(String sequenceType) {
  return "Invalid Sequence Type Error. '$sequenceType' is not a valid sequence type. Please enter the argument 'DNA' or 'RNA' for parameter 'sequenceType'.";
}

String _invalidSequenceErrorMessage(
    String nucleotide, int index, String sequenceType) {
  return "Invalid ${sequenceType.toUpperCase()} Sequence Error. Character '$nucleotide' found at index position $index (zero-based) is not a valid ${sequenceType.toUpperCase()} nucleotide.";
}

Map<String, int> countNucleotideFrequency(
    String nucleotideSequence, sequenceType) {
  Map<String, String> validationMap = validateNucleotideSequence(
      nucleotideSequence = nucleotideSequence, sequenceType = sequenceType);

  String validNucleotideSequence = validationMap[kSequence];

  Map<String, int> nucleotideFrequencyMap = {};
  validNucleotideSequence.split('').forEach((nucleotide) {
    if (nucleotideFrequencyMap.containsKey(nucleotide)) {
      nucleotideFrequencyMap[nucleotide]++;
    } else {
      nucleotideFrequencyMap[nucleotide] = 1;
    }
  });

  return nucleotideFrequencyMap;
}

String _validateSequenceType(String sequenceType) {
  String lowerSequenceType = sequenceType.toLowerCase();

  if (lowerSequenceType != kDNA && lowerSequenceType != kRNA) {
    throw (_invalidSequenceTypeErrorMessage(sequenceType = sequenceType));
  }
  return lowerSequenceType;
}

String generateRandomNucleotideSequence(
    int sequenceLength, String sequenceType) {
  String validSequenceType = _validateSequenceType(sequenceType);

  Random _random = Random();

  if (validSequenceType == kDNA) {
    String DNANucleotidesString = DNANucleotides.join();
    return String.fromCharCodes(
      Iterable.generate(
        sequenceLength,
        (_) => DNANucleotidesString.codeUnitAt(
          _random.nextInt(DNANucleotidesString.length),
        ),
      ),
    );
  }

  String RNANucleotidesString = RNANucleotides.join();
  return String.fromCharCodes(
    Iterable.generate(
      sequenceLength,
      (_) => RNANucleotidesString.codeUnitAt(
        _random.nextInt(RNANucleotidesString.length),
      ),
    ),
  );
}

String transcribeDNAToRNA(String nucleotideSequence) {
  // Fix sequenceType = 'DNA'
  Map<String, String> validationMap =
      validateNucleotideSequence(nucleotideSequence = nucleotideSequence, kDNA);

  String validDNASequence = validationMap[kSequence];

  return validDNASequence.replaceAll(RegExp(r'T'), 'U');
}

String generateComplementaryStrand(
    String nucleotideSequence, String sequenceType, bool reverse) {
  Map<String, String> validationMap = validateNucleotideSequence(
      nucleotideSequence = nucleotideSequence, sequenceType = sequenceType);

  String validNucleotideSequence = validationMap[kSequence];
  String validSequenceType = validationMap[kType];

  String complementarySequence = validNucleotideSequence
      .split('')
      .map((nucleotide) => validSequenceType == kDNA
          ? DNAComplementaryNucleotides[nucleotide]
          : RNAComplementaryNucleotides[nucleotide])
      .join();

  return reverse
      ? _reverseNucleotideSequence(complementarySequence)
      : complementarySequence;
}

String _reverseNucleotideSequence(String nucleotideSequence) {
  return nucleotideSequence.split('').reversed.join('');
}

double calculateGCContent(String nucleotideSequence, String sequenceType) {
  Map<String, String> validationMap = validateNucleotideSequence(
      nucleotideSequence = nucleotideSequence, sequenceType = sequenceType);

  String validNucleotideSequence = validationMap[kSequence];

  int GCCount = 0;
  validNucleotideSequence.split('').forEach((nucleotide) {
    nucleotide == 'G' || nucleotide == 'C' ? GCCount++ : null;
  });

  return num.parse(
      ((GCCount / validNucleotideSequence.length) * 100).toStringAsFixed(2));
}

// TODO - Implement
void detectSequenceType(String nucleotideSequence) {}

int calculateHammingDistance() {}

void main() {
  String DNASequence = generateRandomNucleotideSequence(40, 'DNA');
  String RNASequence = generateRandomNucleotideSequence(40, 'RNA');
  print(validateNucleotideSequence(DNASequence, 'DNA'));
  print(countNucleotideFrequency(DNASequence, 'DNA'));
  print(transcribeDNAToRNA('GATGGAACTTGACTACGTAAATT'));
  print(generateComplementaryStrand('AAAACCCGGT', 'DNA', true));
  print(calculateGCContent('CCACCCTCGTGGTAGGCAGTAGGTGGAAT', 'DNA'));
}
