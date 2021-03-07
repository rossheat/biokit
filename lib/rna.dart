import 'package:meta/meta.dart';
import 'package:biokit/nucleotides.dart';
import 'package:biokit/structs.dart';
import 'dart:math';

class RNA extends Nucleotides {
  RNA({required String seq}) : super(seq: seq, type: 'rna');

  /// Transcribes RNA back to DNA.
  String revTranscribe() => this.seq.replaceAll('U', 'T');

  static RNA random({required int len}) {
    Random _rand = Random();
    String rnaNucsStr = Structs.rnaNucs.join();
    String seq = String.fromCharCodes(
      Iterable.generate(
        len,
        (_) => rnaNucsStr.codeUnitAt(
          _rand.nextInt(rnaNucsStr.length),
        ),
      ),
    );
    return RNA(seq: seq);
  }
}
