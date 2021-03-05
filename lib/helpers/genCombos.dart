import 'package:meta/meta.dart';

List<String> genCombos({@required String seq, sorted = false}) {
  List<String> listSeq = seq.split("");
  List<String> combinations = [];
  for (int i = 0; i < listSeq.length; i++) {
    if (i != listSeq.length - 1) {
      combinations.add(listSeq[i]);
    }
    List<String> temp = [listSeq[i]];
    for (int j = i + 1; j < listSeq.length; j++) {
      temp.add(listSeq[j]);
      combinations.add(temp.join());
    }
  }

  if (sorted) {
    // Sort with longest combination first
    combinations.sort((b, a) => a.length.compareTo(b.length));
    return combinations;
  }
  return combinations;
}
