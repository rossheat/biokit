import 'package:biokit/sequence.dart';
import 'package:biokit/structs.dart';
import 'package:meta/meta.dart';
import 'dart:math';

class Peptide extends Sequence {
  Peptide({@required String seq}) : super(seq: seq, type: 'pep');

  /// Calculates the frequency of each amino acid.
  Map<String, int> freq() => super.freq();

  /// Calculates the Monoisotopic Mass.
  double monoMass({int roundTo = 3}) {
    double totalMonoMass = 0;

    for (var aa in this.seq.split('')) {
      totalMonoMass += Structs.aaToMonoMass[aa];
    }
    return double.parse(totalMonoMass.toStringAsFixed(roundTo));
  }

  static Peptide random({@required int len}) {
    Random _rand = Random();
    String aminoAcidsStr = Structs.aminoAcids.join();
    String seq = String.fromCharCodes(
      Iterable.generate(
        len,
        (_) => aminoAcidsStr.codeUnitAt(
          _rand.nextInt(aminoAcidsStr.length),
        ),
      ),
    );
    return Peptide(seq: seq);
  }
}
