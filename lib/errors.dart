import 'package:meta/meta.dart';

class Errors {
  static String invalidSeq({required String mon, required int idx, required String type}) {
    return "Invalid ${type.toUpperCase()} Sequence Error. Character '$mon' found at index position $idx (zero-based) is not a valid ${type.toUpperCase()} monomer.";
  }

  static String invalidType({required String type}) {
    return "Invalid Sequence Type Error. '$type' is not a valid sequence type. Please enter the argument 'dna', 'rna' or 'pep' for [type].";
  }
}
