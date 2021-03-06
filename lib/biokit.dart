library biokit;

import 'package:biokit/seq.dart';
export 'package:biokit/seq.dart';

void main() async {
  DNA dna = DNA(seq: 'ATCG');
  print(dna.gcContent());
}
