import 'dart:io';
import 'dart:convert';
import 'package:meta/meta.dart';

Future<String> uniprotIdToFASTA({@required String uniprotId}) async {
  Uri uri = Uri.parse('http://www.uniprot.org/uniprot/$uniprotId.fasta');

  var request = await HttpClient().getUrl(uri);

  var response = await request.close();

  await for (var contents in response.transform(Utf8Decoder())) {
    return contents;
  }
  return 'Error retrieving protein with uniprot ID $uniprotId';
}
