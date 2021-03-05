import 'package:meta/meta.dart';

String spliceIntrons({@required nucSeq, @required List<String> introns}) {
  introns.sort((b, a) => a.length.compareTo(b.length));
  for (String intron in introns) {
    nucSeq = nucSeq.replaceAll(intron, '');
  }
  return nucSeq;
}
