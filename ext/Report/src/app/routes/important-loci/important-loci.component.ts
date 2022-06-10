import {Component, Input, OnInit} from '@angular/core';
import {TfDataGetterService} from "../../services/tf-data-getter.service";

@Component({
  selector: 'app-important-loci',
  templateUrl: './important-loci.component.html',
  styleUrls: ['./important-loci.component.css']
})
export class ImportantLociComponent implements OnInit {
  data: {};

  constructor(private dataGetter: TfDataGetterService) {
    this.data = dataGetter.getData().importantLoci;
  }

  ngOnInit(): void {
  }
}
