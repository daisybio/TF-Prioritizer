import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../../types/types";
import {geneCardUrl} from "../../../../assets/constants";

@Component({
  selector: 'cross-links',
  templateUrl: './cross-links.component.html',
  styleUrls: ['./cross-links.component.css']
})
export class CrossLinksComponent implements OnInit {

  @Input()
    // @ts-ignore
  tfGroup: TranscriptionFactorGroup;

  @Input()
  exclude: string = "";

  geneCardUrl: string;

  geneCardDropDownVisible: boolean = false;

  constructor() {
    this.geneCardUrl = geneCardUrl;
  }

  ngOnInit(): void {
  }

  openUrl(url: string) {
    window.open(url);
  }
}
