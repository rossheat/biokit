library biokit;

import 'package:biokit/exports.dart';
export 'package:biokit/exports.dart';

void main() {
  DNA dna = DNA(seq: 'ATCG');
  print(dna.freq());
}
