import 'package:meta/meta.dart';
import '../constants/strings.dart';
import 'validateNucSeq.dart';

Map<String, String> detectSeqType({@required String nucSeq}) {
  bool isValidDNASeq;
  bool isValidRNASeq;
  String invalidDNAErrMsg;
  String invalidRNAErrMsg;

  try {
    validateNucSeq(nucSeq: nucSeq, seqType: kDNA);
    isValidDNASeq = true;
  } catch (errorMessage) {
    isValidDNASeq = false;
    invalidDNAErrMsg = errorMessage;
  }

  try {
    validateNucSeq(nucSeq: nucSeq, seqType: kRNA);
    isValidRNASeq = true;
  } catch (errMsg) {
    isValidRNASeq = false;
    invalidRNAErrMsg = errMsg;
  }

  if (isValidDNASeq == false && isValidRNASeq == false) {
    // Sequence is not DNA or RNA
    return {
      kResult: kInvalidSeq,
      kInvalidDNAErrMsg: invalidDNAErrMsg,
      kInvalidRNAErrMsg: invalidRNAErrMsg
    };
  } else if (isValidDNASeq == true && isValidRNASeq == false) {
    // Sequence is DNA, not RNA
    return {kResult: kDNA, kInvalidRNAErrMsg: invalidRNAErrMsg};
  } else if (isValidDNASeq == false && isValidRNASeq == true) {
    // Sequence is RNA, not DNA
    return {kResult: kRNA, kInvalidDNAErrMsg: invalidDNAErrMsg};
  } else {
    // Cannot tell if sequence is DNA or RNA (No T or U present in sequence)
    return {kResult: kAmbigSeq};
  }
}
