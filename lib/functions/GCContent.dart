import 'package:meta/meta.dart';
import 'package:biokit/helpers/validateNucSeq.dart';
import '../constants/strings.dart';

double gcContent({@required String nucSeq, @required String seqType}) {
  Map<String, String> validationMap =
      validateNucSeq(nucSeq: nucSeq, seqType: seqType);

  String validNucleotideSequence = validationMap[kSeq];

  int gcCount = 0;
  validNucleotideSequence.split('').forEach((nucleotide) {
    nucleotide == 'G' || nucleotide == 'C' ? gcCount++ : null;
  });

  return num.parse(
    ((gcCount / validNucleotideSequence.length) * 100).toStringAsFixed(2),
  );
}
