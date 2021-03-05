import 'package:biokit/constants/maps.dart';
import 'package:biokit/helpers/validateAASeq.dart';
import 'package:meta/meta.dart';

double monoMass({@required String aaSeq, roundTo = 3}) {
  // String validAASeq = validateAASeq(aaSeq: aaSeq);
  double totalMonoMass = 0;

  for (var aa in aaSeq.split('')) {
    totalMonoMass += aaToMonoMass[aa];
  }
  return double.parse(totalMonoMass.toStringAsFixed(roundTo));
}
