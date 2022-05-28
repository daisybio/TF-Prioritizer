import {Component, Input, OnInit} from '@angular/core';
import {geneCardUrl} from "../../../../assets/constants";
import {TranscriptionFactor} from "../../../types/types";
import {OpenExternalUrlService} from "../../../services/open-external-url.service";

@Component({
  selector: 'data-header',
  templateUrl: './data-header.component.html',
  styleUrls: ['./data-header.component.css']
})
export class DataHeaderComponent implements OnInit {
  @Input()
    // @ts-ignore
  title: string;

  @Input()
  information: string | undefined;

  @Input()
    // @ts-ignore
  transcriptionFactor: TranscriptionFactor;

  informationVisible: boolean = false;
  geneCardUrl = geneCardUrl;

  urlOpener: OpenExternalUrlService;

  constructor(urlOpener: OpenExternalUrlService) {
    this.urlOpener = urlOpener;
  }

  ngOnInit(): void {
  }

  toggleInformationVisibility() {
    this.informationVisible = !this.informationVisible;
  }
}
