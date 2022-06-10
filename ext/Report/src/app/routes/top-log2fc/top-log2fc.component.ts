import {Component, OnInit} from '@angular/core';
import {TfDataGetterService} from "../../services/tf-data-getter.service";

@Component({
  selector: 'app-top-log2fc',
  templateUrl: './top-log2fc.component.html',
  styleUrls: ['./top-log2fc.component.css']
})
export class TopLog2fcComponent implements OnInit {
  data: {};

  constructor(private dataGetter: TfDataGetterService) {
    this.data = dataGetter.getData().topLog2fc;
  }

  ngOnInit(): void {
  }
}
