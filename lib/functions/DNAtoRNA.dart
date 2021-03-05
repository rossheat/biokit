import 'package:meta/meta.dart';
import 'package:biokit/helpers/validateNucSeq.dart';
import '../constants/strings.dart';

String dnaToRNA({@required String dnaSeq}) {
  Map<String, String> validationMap =
      validateNucSeq(nucSeq: dnaSeq, seqType: kDNA);
  String validDNASeq = validationMap[kSeq];
  return validDNASeq.replaceAll(RegExp(r'T'), 'U');
}
