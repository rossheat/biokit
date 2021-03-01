import 'dart:math';
import '../constants/lists.dart';
import '../constants/strings.dart';
import 'validateSequenceType.dart';

String makeNucleotideSequence(int sequenceLength, String sequenceType) {
  String validSequenceType = validateSequenceType(sequenceType);

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
