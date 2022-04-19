import {Component, Input, OnInit} from '@angular/core';

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

  informationVisible: boolean = false;

  constructor() {
  }

  ngOnInit(): void {
  }

  toggleInformationVisibility() {
    this.informationVisible = !this.informationVisible;
  }
}
