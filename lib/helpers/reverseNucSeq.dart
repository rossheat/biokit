import 'package:meta/meta.dart';

String reverseNucSeq({@required String nucSeq}) {
  return nucSeq.split('').reversed.join('');
}
