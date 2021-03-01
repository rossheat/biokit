import '../constants/strings.dart';
import '../errors/invalidSequenceTypeErrorMessage.dart';

String validateSequenceType(String sequenceType) {
  String lowerSequenceType = sequenceType.toLowerCase();

  if (lowerSequenceType != kDNA && lowerSequenceType != kRNA) {
    throw (invalidSequenceTypeErrorMessage(sequenceType = sequenceType));
  }
  return lowerSequenceType;
}
