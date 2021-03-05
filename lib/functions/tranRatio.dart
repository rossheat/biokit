import 'package:biokit/constants/lists.dart';
import 'package:meta/meta.dart';

double tranRatio({@required String dnaSeqOne, @required String dnaSeqTwo}) {
  //
  if (dnaSeqOne.length != dnaSeqTwo.length) {
    throw ('Unequal Sequence Lengths Error. ');
  }

  int transitionCount = 0;
  int transversionCount = 0;

  dnaSeqOne.split('').asMap().forEach((index, seqOneNuc) {
    if (seqOneNuc != dnaSeqTwo[index]) {
      if (dnaTransitions.contains(seqOneNuc + dnaSeqTwo[index])) {
        transitionCount++;
      } else if (dnaTransversions.contains(seqOneNuc + dnaSeqTwo[index])) {
        transversionCount++;
      }
    }
  });

  return transitionCount / transversionCount;
}
