import '../constants/strings.dart';
import '../helpers/validateNucleotideSequence.dart';

double GCContent(String nucleotideSequence, String sequenceType) {
  Map<String, String> validationMap = validateNucleotideSequence(
      nucleotideSequence = nucleotideSequence, sequenceType = sequenceType);

  String validNucleotideSequence = validationMap[kSequence];

  int GCCount = 0;
  validNucleotideSequence.split('').forEach((nucleotide) {
    nucleotide == 'G' || nucleotide == 'C' ? GCCount++ : null;
  });

  return num.parse(
    ((GCCount / validNucleotideSequence.length) * 100).toStringAsFixed(2),
  );
}
